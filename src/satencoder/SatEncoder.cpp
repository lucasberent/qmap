#include "satencoder/SatEncoder.hpp"

#include "CircuitOptimizer.hpp"
#include "boost/dynamic_bitset.hpp"
#include "exact/ExactMapper.hpp"
#include "z3++.h"

#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.
#include <chrono>
//
// Created by lucas on 25/01/2022.
//
bool SatEncoder::testEqual(qc::QuantumComputation& circuitOne, qc::QuantumComputation& circuitTwo, std::vector<std::string>& inputs) {
    if (!isClifford(circuitOne) || !isClifford(circuitTwo)) {
        std::cerr << "Circuits are not Clifford circuits" << std::endl;
        return false;
    }
    qc::DAG                           dagOne     = qc::CircuitOptimizer::constructDAG(circuitOne);
    qc::DAG                           dagTwo     = qc::CircuitOptimizer::constructDAG(circuitTwo);
    SatEncoder::CircuitRepresentation circOneRep = preprocessCircuit(dagOne, inputs);
    SatEncoder::CircuitRepresentation circTwoRep = preprocessCircuit(dagTwo, inputs);

    const auto  before = std::chrono::high_resolution_clock::now();
    z3::context ctx{};
    z3::solver  solver(ctx);
    constructMiterInstance(circOneRep, circTwoRep, solver);
    const auto after                   = std::chrono::high_resolution_clock::now();
    const auto satConstructionDuration = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
    std::cout << "SAT construction complete - elapsed time (ms) for this task: " << satConstructionDuration.count() << std::endl;
    return !isSatisfiable(solver); // if unsat circuits are equal
    //todo print stats
}
void SatEncoder::checkSatisfiability(qc::QuantumComputation& circuitOne, std::vector<std::string>& inputs) {
    if (!isClifford(circuitOne)) {
        std::cerr << "Circuit is not Clifford Circuit." << std::endl;
        return;
    }
    qc::DAG    dag                   = qc::CircuitOptimizer::constructDAG(circuitOne);
    auto       before                = std::chrono::high_resolution_clock::now();
    auto       circRep               = preprocessCircuit(dag, inputs);
    auto       after                 = std::chrono::high_resolution_clock::now();
    const auto preprocessingDuration = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
    std::cout << "Preprocessing construction complete - elapsed time (ms) for this task: " << preprocessingDuration.count() << std::endl;

    before = std::chrono::high_resolution_clock::now();
    z3::context ctx{};
    z3::solver  solver(ctx);
    constructSatInstance(circRep, solver);
    after                        = std::chrono::high_resolution_clock::now();
    auto satConstructionDuration = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
    std::cout << "SAT construction complete - elapsed time (ms) for this task: " << satConstructionDuration.count() << std::endl;

    // todo print to file
    // data: ckt size: #gates/depth, #different generators (states), #sat variables, time needed for construction
    // ckt name, ckt size, equiv of not, #sat variables, #different generators
    // in z3 stats: :mk-bool-var :mk-binary-clause :mk-ternary-clause :mk-clause
    //print stats;
    isSatisfiable(solver);
}

bool SatEncoder::isSatisfiable(solver& solver) {
    bool result            = false;
    auto before            = std::chrono::high_resolution_clock::now();
    auto sat               = solver.check();
    auto after             = std::chrono::high_resolution_clock::now();
    auto z3SolvingDuration = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
    std::cout << "Z3 solving complete - elapsed time (ms) for this task: " << z3SolvingDuration.count() << std::endl;
    if (sat == check_result::sat) {
        std::cout << "SATISFIABLE" << std::endl;
        result = true;
    } else if (sat == check_result::unsat) {
        std::cout << "UNSATISFIABLE" << std::endl;
    } else {
        std::cerr << "UNKNOWN" << std::endl;
    }
    std::cout << "STATS: " << solver.statistics();
    //std::cout << "STATS: " << solver.get_model();
    return result;
}

SatEncoder::CircuitRepresentation SatEncoder::preprocessCircuit(qc::DAG& dag, std::vector<std::string>& inputs) {
    unsigned long                     inputSize      = dag.size();
    unsigned long                     nrOfLevels     = 0;
    std::size_t                       nrOfOpsOnQubit = 0;
    unsigned long                     tmp;
    unsigned long                     uniqueGenCnt = 0;
    std::vector<QState>               states;
    SatEncoder::CircuitRepresentation representation;

    unsigned long nrOfQubits = dag.size();
    //compute nr of levels of ckt = #generators needed per input state
    for (std::size_t i = 0U; i < inputSize; i++) {
        tmp = dag.at(i).size();
        if (tmp > nrOfLevels) {
            nrOfLevels = tmp;
        }
    }
    std::vector<std::map<boost::uuids::uuid, boost::uuids::uuid>> vec(nrOfLevels);
    representation.generatorMappings = vec;

    if (inputs.size() > 1) {
        for (auto& input: inputs) {
            states.push_back(initializeState(nrOfQubits, input));
        }
    } else {
        states.push_back(initializeState(nrOfQubits, {}));
    }

    // store generators of input state
    for (auto& state: states) {
        auto               initLevelGenerator = state.getLevelGenerator();
        auto               inspair            = generators.emplace(initLevelGenerator, uniqueGenCnt); // put generator into global map if not already present
        boost::uuids::uuid genId{};
        if (inspair.second) { // if a new generator has been computed by this level (i.e., state changed)
            uniqueGenCnt++;
            genId = boost::uuids::random_generator()();
            generatorIdMap.emplace(initLevelGenerator, genId);
        } else {
            genId = generatorIdMap.at(initLevelGenerator);
        }
        representation.idGeneratorMap.emplace(genId, initLevelGenerator);
        state.SetPrevGenId(genId);
        //no generator <> generator mapping for initial level
        std::cout << "Init State:" << std::endl;
        state.printStateTableau();
    }

    for (std::size_t levelCnt = 0; levelCnt < nrOfLevels; levelCnt++) {
        for (std::size_t qubitCnt = 0U; qubitCnt < inputSize; qubitCnt++) { //apply operation of current level for each qubit
            nrOfOpsOnQubit = dag.at(qubitCnt).size();

            if (levelCnt < nrOfOpsOnQubit) {
                if (!dag.at(qubitCnt).empty() && dag.at(qubitCnt).at(levelCnt) != nullptr) {
                    stats.nrOfGates++;
                    auto          gate    = dag.at(qubitCnt).at(levelCnt)->get();
                    unsigned long target  = gate->getTargets().at(0U);          //we assume we only have 1 target
                    unsigned long control = gate->getControls().begin()->qubit; //we assume we only have 1 control

                    //apply gate of level to each generator
                    for (auto& currState: states) {
                        if (gate->getType() == qc::OpType::H) {
                            std::cout << "level " << levelCnt << " applying H to qubit " << target << std::endl;
                            currState.applyH(target);
                        } else if (gate->getType() == qc::OpType::S) {
                            std::cout << "level " << levelCnt << " applying S to qubit " << target << std::endl;
                            currState.applyS(target);
                        } else if (gate->isControlled() && gate->getType() == qc::OpType::X) { //CNOT
                            if (qubitCnt == control) {                                         //CNOT is for control and target in DAG, only apply if current qubit is control
                                std::cout << "level " << levelCnt << " applying CNOT with target: " << target << " control: " << control << std::endl;
                                currState.applyCNOT(control, target);
                            }
                        } else {
                            //todo error operation not supported
                        }
                    }
                }
            }
        }
        for (auto& state: states) {
            std::cout << "state " << std::endl;
            state.printStateTableau();
            auto               currLevelGen = state.getLevelGenerator();                      //extract generator representation from tableau (= list of paulis for each qubit)
            auto               inspair      = generators.emplace(currLevelGen, uniqueGenCnt); // put generator into global map if not already present
            boost::uuids::uuid genId{};
            if (inspair.second) { // new generator because newly inserted
                uniqueGenCnt++;
                genId = boost::uuids::random_generator()();
                generatorIdMap.emplace(currLevelGen, genId); // global generator <> id mapping for quick reverse lookup
            } else {                                         // generator already in global map
                genId = generatorIdMap.at(currLevelGen);
            }
            representation.idGeneratorMap.emplace(genId, currLevelGen);                                        // id <-> generator mapping
            representation.generatorMappings.at(levelCnt).insert(std::make_pair(state.GetPrevGenId(), genId)); // generator <> generator mapping at position level in list
            state.SetPrevGenId(genId);                                                                         // update previous generator id for next state
        }
    }
    return representation;
}
//construct z3 instance from preprocessing information
void SatEncoder::constructSatInstance(SatEncoder::CircuitRepresentation& circuitRepresentation, z3::solver& solver) {
    // number of unique generators that need to be encoded
    const auto generatorCnt = generators.size();
    stats.nrOfGenerators    = generatorCnt;

    // bitwidth required to encode the generators
    const auto bitwidth = static_cast<std::size_t>(std::ceil(std::log2(generatorCnt)));

    // whether the number of generators is a power of two or not
    bool blockingConstraintsNeeded = std::log2(generatorCnt) < static_cast<double>(bitwidth);

    // z3 context used throughout this function
    auto& ctx = solver.ctx();

    const auto        depth = circuitRepresentation.generatorMappings.size();
    std::vector<expr> vars{};
    vars.reserve(depth + 1U);
    std::string bvName = "x^";

    std::cout << "Variables: " << std::endl;
    for (std::size_t k = 0U; k <= depth; k++) {
        // create bitvector [x^k]_2 with respective bitwidth for each level k of ckt
        std::stringstream ss{};
        ss << bvName << k; //
        vars.emplace_back(ctx.bv_const(ss.str().c_str(), bitwidth));
        std::cout << vars.back() << std::endl;
    }

    std::cout << "Functional Constraints: " << std::endl;
    for (std::size_t i = 0U; i < depth; i++) {
        const auto layer = circuitRepresentation.generatorMappings.at(i); // generator<>generator map for level i
        for (const auto& [from, to]: layer) {
            const auto g1 = generators.at(circuitRepresentation.idGeneratorMap.at(from));
            const auto g2 = generators.at(circuitRepresentation.idGeneratorMap.at(to));

            // create [x^l]_2 = i => [x^l']_2 = k for each generator mapping
            const auto left  = vars[i] == ctx.bv_val(static_cast<std::uint64_t>(g1), bitwidth);
            const auto right = vars[i + 1U] == ctx.bv_val(static_cast<std::uint64_t>(g2), bitwidth);
            const auto cons  = implies(left, right);
            solver.add(cons);

            stats.nrOfFunctionalConstr++;
            std::cout << cons << std::endl;
        }
    }

    if (blockingConstraintsNeeded) {
        std::cout << "Blocking Constraints: " << std::endl;
        for (const auto& var: vars) {
            const auto cons = var < ctx.bv_val(static_cast<std::uint64_t>(generatorCnt), bitwidth);
            solver.add(cons); // [x^l]_2 < m
            std::cout << cons << std::endl;
        }
    }
}

void SatEncoder::constructMiterInstance(SatEncoder::CircuitRepresentation& circOneRep, SatEncoder::CircuitRepresentation& circTwoRep, z3::solver& solver) {
    // number of unique generators that need to be encoded
    const auto generatorCnt = generators.size();
    stats.nrOfGenerators    = generatorCnt;

    // bitwidth required to encode the generators
    const auto bitwidth = static_cast<std::size_t>(std::ceil(std::log2(generatorCnt)));

    // whether the number of generators is a power of two or not
    bool blockingConstraintsNeeded = std::log2(generatorCnt) < static_cast<double>(bitwidth);

    // z3 context used throughout this function
    auto& ctx = solver.ctx();

    /// encode first circuit
    const auto        depthOne = circOneRep.generatorMappings.size();
    std::vector<expr> varsOne{};
    varsOne.reserve(depthOne + 1U);
    std::string bvName = "x^";

    std::cout << "Variables: " << std::endl;
    for (std::size_t k = 0U; k <= depthOne; k++) {
        // create bitvector [x^k]_2 with respective bitwidth for each level k of ckt
        std::stringstream ss{};
        ss << bvName << k; //
        varsOne.emplace_back(ctx.bv_const(ss.str().c_str(), bitwidth));
        std::cout << varsOne.back() << std::endl;
    }

    std::cout << "Functional Constraints: " << std::endl;
    for (std::size_t i = 0U; i < depthOne; i++) {
        const auto layer = circOneRep.generatorMappings.at(i); // generator<>generator map for level i
        for (const auto& [from, to]: layer) {
            const auto g1 = generators.at(circOneRep.idGeneratorMap.at(from));
            const auto g2 = generators.at(circOneRep.idGeneratorMap.at(to));

            // create [x^l]_2 = i <=> [x^l']_2 = k for each generator mapping
            const auto left    = varsOne[i] == ctx.bv_val(static_cast<std::uint64_t>(g1), bitwidth);
            const auto right   = varsOne[i + 1U] == ctx.bv_val(static_cast<std::uint64_t>(g2), bitwidth);
            const auto cons    = implies(left, right);
            const auto revcons = implies(right, left);
            solver.add(cons);
            solver.add(revcons);
            stats.nrOfFunctionalConstr++;
            std::cout << cons << std::endl;
            std::cout << revcons << std::endl;
        }
    }

    if (blockingConstraintsNeeded) {
        std::cout << "Blocking Constraints: " << std::endl;
        for (const auto& var: varsOne) {
            const auto cons = var < ctx.bv_val(static_cast<std::uint64_t>(generatorCnt), bitwidth);
            solver.add(cons); // [x^l]_2 < m
            std::cout << cons << std::endl;
        }
    }

    /// encode second circuit
    auto              depthTwo = circTwoRep.generatorMappings.size();
    std::vector<expr> varsTwo{};
    varsOne.reserve(depthTwo + 1U);
    bvName = "x'^";

    std::cout << "Variables: " << std::endl;
    for (std::size_t k = 0U; k <= depthTwo; k++) {
        // create bitvector [x^k]_2 with respective bitwidth for each level k of ckt
        std::stringstream ss{};
        ss << bvName << k; //
        varsTwo.emplace_back(ctx.bv_const(ss.str().c_str(), bitwidth));
        std::cout << varsTwo.back() << std::endl;
    }

    std::cout << "Functional Constraints: " << std::endl;
    for (std::size_t i = 0U; i < depthTwo; i++) {
        const auto layer = circTwoRep.generatorMappings.at(i); // generator<>generator map for level i
        for (const auto& [from, to]: layer) {
            const auto g1 = generators.at(circTwoRep.idGeneratorMap.at(from));
            const auto g2 = generators.at(circTwoRep.idGeneratorMap.at(to));

            // create [x^l]_2 = i <=> [x^l']_2 = k for each generator mapping
            const auto left    = varsTwo[i] == ctx.bv_val(static_cast<std::uint64_t>(g1), bitwidth);
            const auto right   = varsTwo[i + 1U] == ctx.bv_val(static_cast<std::uint64_t>(g2), bitwidth);
            const auto cons    = implies(left, right);
            const auto revcons = implies(right, left);
            solver.add(revcons);
            solver.add(cons);

            stats.nrOfFunctionalConstr++;
            std::cout << cons << std::endl;
            std::cout << revcons << std::endl;
        }
    }

    if (blockingConstraintsNeeded) {
        std::cout << "Blocking Constraints: " << std::endl;
        for (const auto& var: varsTwo) {
            const auto cons = var < ctx.bv_val(static_cast<std::uint64_t>(generatorCnt), bitwidth);
            solver.add(cons); // [x^l]_2 < m
            std::cout << cons << std::endl;
        }
    }

    // create miter structure
    // if initial signals are the same, then the final signals have to be equal as well
    const auto initial = varsOne.front() == varsTwo.front();
    const auto final   = varsOne.back() != varsTwo.back();
    const auto miter   = implies(initial, final);
    std::cout << "Miter: " << miter << std::endl;
    solver.add(miter);
}

bool SatEncoder::isClifford(qc::QuantumComputation& qc) {
    qc::OpType opType;
    for (const auto& op: qc) {
        opType = op->getType();
        if (opType != qc::OpType::H &&
            opType != qc::OpType::S &&
            opType != qc::OpType::X &&
            opType != qc::OpType::Z &&
            opType != qc::OpType::Y &&
            opType != qc::OpType::I) {
            return false;
        }
    }
    return true;
}
std::vector<boost::dynamic_bitset<>> SatEncoder::QState::getLevelGenerator() const {
    std::size_t                          size = (2U * n) + 1U;
    std::vector<boost::dynamic_bitset<>> result{};

    for (std::size_t i = 0U; i < n; i++) {
        boost::dynamic_bitset<> gen(size);
        for (std::size_t j = 0; j < n; j++) {
            gen[j] = x.at(i)[j];
        }
        for (std::size_t j = 0; j < n; j++) {
            gen[n + j] = z.at(i)[j];
        }
        if (r.at(i) == 1) {
            gen[n + n] = true;
        } else {
            gen[n + n] = false; //either 0 or 1 possible for phase
        }
        result.emplace_back(gen);
    }

    return result;
}
SatEncoder::QState SatEncoder::initializeState(unsigned long nrOfQubits, std::string input) {
    SatEncoder::QState result;
    result.SetN(nrOfQubits);
    std::vector<boost::dynamic_bitset<>> tx(nrOfQubits);
    std::vector<boost::dynamic_bitset<>> tz(nrOfQubits);
    std::vector<int>                     tr(nrOfQubits, 0);

    for (std::size_t i = 0U; i < nrOfQubits; i++) {
        tx[i] = boost::dynamic_bitset(nrOfQubits);
        tz[i] = boost::dynamic_bitset(nrOfQubits);
        for (std::size_t j = 0; j < nrOfQubits; j++) {
            if (i == j) {
                tz[i][j] = true; // initial 0..0 state corresponds to x matrix all zero and z matrix = Id_n
            }
        }
    }

    result.SetX(tx);
    result.SetZ(tz);
    result.SetR(tr);
    if (!input.empty()) { //
        for (std::size_t i = 0U; i < input.length(); i++) {
            switch (input[i]) {
                case 'Z': // stab by -Z = |1>
                    result.applyH(i);
                    result.applyS(i);
                    result.applyS(i);
                    result.applyH(i);
                    break;
                case 'x': // stab by X = |+>
                    result.applyH(i);
                    break;
                case 'X': // stab by -X = |->
                    result.applyH(i);
                    result.applyS(i);
                    result.applyS(i);
                    break;
                case 'y': // stab by Y = |0> + i|1>
                    result.applyH(i);
                    result.applyS(i);
                    break;
                case 'Y': // stab by -Y = |0> - i|1>
                    result.applyH(i);
                    result.applyS(i);
                    result.applyS(i);
                    result.applyS(i);
            }
        }
    }
    return result;
}

void SatEncoder::QState::applyCNOT(unsigned long control, unsigned long target) {
    if (target > n || target < 0 ||
        control > n || control < 0) {
        return;
    }
    for (std::size_t i = 0U; i < n; ++i) {
        r[i] ^= (x[i][control] * z[i][target]) * (x[i][target] ^ z[i][control] ^ 1);
        x[i][target] ^= x[i][control];
        z[i][control] ^= z[i][target];
    }
}
void SatEncoder::QState::applyH(unsigned long target) {
    if (target > n || target < 0) {
        return;
    }
    for (std::size_t i = 0U; i < n; i++) {
        r[i] ^= x[i][target] * z[i][target];
        x[i][target] ^= z[i][target];
        z[i][target] = x[i][target] ^ z[i][target];
        x[i][target] ^= z[i][target];
    }
}
void SatEncoder::QState::applyS(unsigned long target) {
    if (target > n || target < 0) {
        return;
    }
    for (std::size_t i = 0U; i < n; ++i) {
        r[i] ^= x[i][target] * z[i][target];
        z[i][target] ^= x[i][target];
    }
}
void SatEncoder::QState::SetN(unsigned long N) {
    QState::n = N;
}
const std::vector<boost::dynamic_bitset<>>& SatEncoder::QState::GetX() const {
    return x;
}
void SatEncoder::QState::SetX(const std::vector<boost::dynamic_bitset<>>& X) {
    QState::x = X;
}
void SatEncoder::QState::SetZ(const std::vector<boost::dynamic_bitset<>>& Z) {
    QState::z = Z;
}
void SatEncoder::QState::SetR(const std::vector<int>& R) {
    QState::r = R;
}
const boost::uuids::uuid& SatEncoder::QState::GetPrevGenId() const {
    return prevGenId;
}
void SatEncoder::QState::SetPrevGenId(const boost::uuids::uuid& prev_gen_id) {
    prevGenId = prev_gen_id;
}
void SatEncoder::QState::printStateTableau() {
    std::cout << std::endl;
    for (std::size_t i = 0U; i < n; i++) {
        for (std::size_t j = 0; j < n; j++) {
            std::cout << x.at(i)[j];
        }
        std::cout << "|";
        for (std::size_t j = 0; j < n; j++) {
            std::cout << z.at(i)[j];
        }
        std::cout << "|";
        std::cout << r.at(i);
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
