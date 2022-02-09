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
        std::cerr << "Circits are not clifford circuits" << std::endl;
        return false;
    }
    qc::DAG                           dagOne                  = qc::CircuitOptimizer::constructDAG(circuitOne);
    qc::DAG                           dagTwo                  = qc::CircuitOptimizer::constructDAG(circuitTwo);
    SatEncoder::CircuitRepresentation circOneRep              = preprocessCircuit(dagOne, inputs);
    SatEncoder::CircuitRepresentation circTwoRep              = preprocessCircuit(dagTwo, inputs);
    auto                              before                  = std::chrono::high_resolution_clock::now();
    auto                              solver                  = constructMiterInstance(circOneRep, circTwoRep);
    auto                              after                   = std::chrono::high_resolution_clock::now();
    auto                              satConstructionDuration = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
    std::cout << "SAT construction complete - elapsed time for this task: " << satConstructionDuration.count();
    return runZ3(solver);
    //todo print stats
}
void SatEncoder::checkSatisfiability(qc::QuantumComputation& circuitOne, std::vector<std::string>& inputs) {
    if (!isClifford(circuitOne)) {
        std::cerr << "Circuit is not Clifford Circuit." << std::endl;
        return;
    }
    qc::DAG dag                   = qc::CircuitOptimizer::constructDAG(circuitOne);
    auto    before                = std::chrono::high_resolution_clock::now();
    auto    circRep               = this->preprocessCircuit(dag, inputs);
    auto    after                 = std::chrono::high_resolution_clock::now();
    auto    preprocessingDuration = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
    std::cout << "Preprocessing construction complete - elapsed time (ms) for this task: " << preprocessingDuration.count() << std::endl;
    before                       = std::chrono::high_resolution_clock::now();
    auto solver                  = this->constructSatInstance(circRep);
    after                        = std::chrono::high_resolution_clock::now();
    auto satConstructionDuration = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
    std::cout << "SAT construction complete - elapsed time (ms) for this task: " << satConstructionDuration.count() << std::endl;
    // todo print to file
    // data: ckt size: #gates/depth, #different generators (states), #sat variables, time needed for construction
    // ckt name, ckt size, equiv of not, #sat variables, #different generators
    // in z3 stats: :mk-bool-var :mk-binary-clause :mk-ternary-clause :mk-clause
    //print stats;
    runZ3(solver);
}

bool SatEncoder::runZ3(solver& solver) {
    bool result            = false;
    auto before            = std::chrono::high_resolution_clock::now();
    auto sat               = solver.check();
    auto after             = std::chrono::high_resolution_clock::now();
    auto z3SolvingDuration = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
    std::cout << "Z3 solving complete - elapsed time for this task: " << z3SolvingDuration.count();
    if (sat == check_result::sat) {
        std::cout << "SATISFIABLE" << std::endl;
        result = true;
    } else if (sat == check_result::unsat) {
        std::cout << "UNSATISFIABLE" << std::endl;
    } else {
        std::cerr << "UNKNOWN" << std::endl;
    }
    std::cout << "STATS: " << solver.statistics();
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
        auto initLevelGenerator = state.getLevelGenerator();
        auto inspair            = generators.emplace(initLevelGenerator, uniqueGenCnt); // put generator into global map if not already present
        if (inspair.second) {                                                           // if a new generator has been computed by this level (i.e., state changed)
            uniqueGenCnt++;
            boost::uuids::uuid genId = boost::uuids::random_generator()();
            representation.idGeneratorMap.emplace(genId, initLevelGenerator);
            state.SetPrevGenId(genId);
            //no generator <> generator mapping for initial level
            std::cout << "Init State:" << std::endl;
            state.printStateTableau();
        }
    }

    for (std::size_t levelCnt = 0; levelCnt < nrOfLevels; levelCnt++) {
        for (std::size_t qubitCnt = 0U; qubitCnt < inputSize; qubitCnt++) { //apply operation of current level for each qubit
            nrOfOpsOnQubit = dag.at(qubitCnt).size();

            if (levelCnt < nrOfOpsOnQubit) {
                if (!dag.at(qubitCnt).empty() && dag.at(qubitCnt).at(levelCnt) != nullptr) {
                    stats.nrOfGates++;
                    auto          gate    = dag.at(qubitCnt).at(levelCnt)->get();
                    unsigned long target  = gate->getTargets().at(0U);           //we assume we only have 1 target
                    unsigned long control = gate->getControls().begin()->qubit; //we assume we only have 1 control

                    //apply gate of level to each generator
                    for (auto& currState: states) {
                        if (gate->getType() == qc::OpType::H) {
                            std::cout << "applying H to qubit " << target << std::endl;
                            currState.applyH(target);
                        } else if (gate->getType() == qc::OpType::S) {
                            std::cout << "applying S to qubit " << target << std::endl;
                            currState.applyS(target);
                        } else if (gate->isControlled() && gate->getType() == qc::OpType::X) { //CNOT
                            if (qubitCnt == control) {                                         //CNOT is for control and target in DAG, only apply if current qubit is control
                                std::cout << "applying CNOT with target: " << target << " control: " << control << std::endl;
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
            auto currLevelGen = state.getLevelGenerator();                      //extract generator representation from tableau (= list of paulis for each qubit)
            auto inspair      = generators.emplace(currLevelGen, uniqueGenCnt); // put generator into global map if not already present
            if (inspair.second) {                                               // new generator because newly inserted
                uniqueGenCnt++;
                boost::uuids::uuid genId = boost::uuids::random_generator()();
                representation.idGeneratorMap.emplace(genId, currLevelGen);                                        // id <-> generator mapping
                representation.generatorMappings.at(levelCnt).insert(std::make_pair(state.GetPrevGenId(), genId)); // insert generator <> generator mapping at position level in list
                state.SetPrevGenId(genId);                                                                         // update previous generator id for next state
            }
        }
    }
    return representation;
}
//construct z3 instance from preprocessing information
solver SatEncoder::constructSatInstance(SatEncoder::CircuitRepresentation& circuitRepresentation) {
    auto        generatorCnt              = generators.size();
    unsigned    nrOfVariablesNeeded       = ceil(log(generatorCnt));
    auto        depth                     = circuitRepresentation.generatorMappings.size();
    bool        blockingConstraintsNeeded = (generatorCnt < (1ULL << nrOfVariablesNeeded));
    z3::context z3Context;
    //for stats
    //nrOfGenerators = generatorCnt;
    //nrOfSatVars    = nrOfVariablesNeeded;

    // construct variable encoding for generators
    std::vector<expr> satVariables{};
    satVariables.reserve(depth);
    std::string bitvName = "x^";
    std::cout << "Variables: " << std::endl;
    for (std::size_t k = 0; k < depth; k++) {
        std::stringstream varName{};
        // create bitvector [x^k]_2 of size nrOfVariablesNeeded for each level k of ckt
        varName << bitvName << k; //
        satVariables.emplace_back(z3Context.bv_const(varName.str().c_str(), nrOfVariablesNeeded));
        std::cout << satVariables[k] << std::endl;
    }
    z3::solver z3Solver(z3Context);
    std::cout << "Functional Constraints: " << std::endl;
    for (std::size_t i = 1; i < circuitRepresentation.generatorMappings.size(); i++) {
        auto tmp = circuitRepresentation.generatorMappings.at(i); // generator<>generator map for level i
        for (auto& ti: tmp) {
            const auto g1 = generators.at(circuitRepresentation.idGeneratorMap.at(ti.first));
            const auto g2 = generators.at(circuitRepresentation.idGeneratorMap.at(ti.second));
            stats.nrOfFunctionalConstr++;
            auto cons = implies(satVariables[i - 1] == z3Context.bv_val(static_cast<std::uint64_t>(g1), nrOfVariablesNeeded),
                                satVariables[i] == z3Context.bv_val(static_cast<std::uint64_t>(g2), nrOfVariablesNeeded));
            z3Solver.add(cons); // create [x^l]_2 = i => [x^l']_2 = k for each generator mapping
            std::cout << cons << std::endl;
        }
    }
    std::cout << "Blocking Constraints: " << std::endl;
    if (blockingConstraintsNeeded) {
        // each bitvector assignment x^g must be < generatorCnt
        for (auto& satVariable: satVariables) {
            auto bv = ult(satVariable, static_cast<int>(generatorCnt));
            std::cout << bv << std::endl;
            z3Solver.add(bv); // [x^g]_2 < m
        }
    }
    return z3Solver;
}

solver SatEncoder::constructMiterInstance(SatEncoder::CircuitRepresentation& circOneRep, SatEncoder::CircuitRepresentation& circTwoRep) {
    unsigned long generatorCnt              = generators.size();
    stats.nrOfGenerators                    = generatorCnt;
    unsigned long nrOfVariablesNeeded       = ceil(log(generatorCnt));
    auto          depth                     = circOneRep.generatorMappings.size();
    bool          blockingConstraintsNeeded = (generatorCnt < (1ULL << nrOfVariablesNeeded));
    context       z3Context;
    //for stats
    //nrOfGenerators = generatorCnt;
    //nrOfSatVars    = nrOfVariablesNeeded;

    std::vector<expr> satVariables{};
    satVariables.reserve(depth);
    std::string bitvName = "x^";
    std::cout << "Variables: " << std::endl;

    for (std::size_t k = 0; k < depth; k++) {
        // create bitvector [x^k]_2 of size nrOfVariablesNeeded for each level k of ckt
        std::stringstream varName{};
        varName << bitvName << k; //
        satVariables[k] = (z3Context.bv_const(varName.str().c_str(), nrOfVariablesNeeded));
        std::cout << satVariables[k] << std::endl;
    }
    solver z3Solver(z3Context);
    std::cout << "Functional Constraints: " << std::endl;
    for (std::size_t i = 1; i < circOneRep.generatorMappings.size(); i++) {
        auto tmp = circOneRep.generatorMappings.at(i); // generator<>generator map for level i
        for (auto& ti: tmp) {
            auto g1 = generators.at(circOneRep.idGeneratorMap.at(ti.first));
            auto g2 = generators.at(circOneRep.idGeneratorMap.at(ti.second));
            //nrOfFunctionalConstr++; // for stats
            auto cons = implies(satVariables[i - 1] == z3Context.bv_val(static_cast<std::uint64_t>(g1), nrOfVariablesNeeded),
                                satVariables[i] == z3Context.bv_val(static_cast<std::uint64_t>(g2), nrOfVariablesNeeded));
            z3Solver.add(cons); // create [x^l]_2 = i => [x^l']_2 = k for each generator mapping
            std::cout << cons << std::endl;
        }
    }

    if (blockingConstraintsNeeded) {
        for (auto& satVariable: satVariables) {
            z3Solver.add(ult(satVariable, static_cast<int>(generatorCnt))); // [x^l]_2 < m
        }
    }
    // circuit 2
    depth                     = circTwoRep.generatorMappings.size();
    blockingConstraintsNeeded = (generatorCnt < (1ULL << nrOfVariablesNeeded));
    std::cout << "Functional Constraints Circuit 2: " << std::endl;
    // construct functional encoding for gates
    for (std::size_t i = 1; i < circTwoRep.generatorMappings.size(); i++) {
        auto tmp = circTwoRep.generatorMappings.at(i); // generator<>generator map for level i
        for (auto& ti: tmp) {
            auto g1 = generators.at(circTwoRep.idGeneratorMap.at(ti.first));
            auto g2 = generators.at(circTwoRep.idGeneratorMap.at(ti.second));
            //nrOfFunctionalConstr++; // for stats
            auto cons = implies(satVariables[i - 1] == z3Context.bv_val(static_cast<std::uint64_t>(g1), nrOfVariablesNeeded),
                                satVariables[i] == z3Context.bv_val(static_cast<std::uint64_t>(g2), nrOfVariablesNeeded));
            z3Solver.add(cons); // create [x^l]_2 = i => [x^l']_2 = k for each generator mapping
            std::cout << cons << std::endl;
        }
    }

    // miter: check if [g^l]_2 == [g^l']_2 for [g^l]_2 encoding last level of circut one and [g^l']_2 encoding last level of circuit two
    // that is, if generator of last level corresponding to some variable assignment of the variables of the signal of the last level
    // are equal then the 'output' is equal.
    auto maxDepthCircOne = circOneRep.generatorMappings.size() - 1U;
    auto maxDepthCircTwo = circTwoRep.generatorMappings.size() - 1U;
    auto miter           = ugt(satVariables[maxDepthCircOne] ^ satVariables[maxDepthCircTwo], 1);
    std::cout << "Miter: " << miter << std::endl;
    z3Solver.add(miter);
    return z3Solver;
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
    std::vector<boost::dynamic_bitset<>> result(n);
    std::size_t                          idx = 0U;
    for (std::size_t i = 0U; i < n; i++) {
        boost::dynamic_bitset<> gen(size);

        //copy x and z vectors and r bit into one bitvector = 1 generator
        for (std::size_t j = 0; j < n; j++) {
            if (idx < this->n) { //x
                gen[idx++] = this->GetX().at(i)[j];
            } else if (idx < 2 * this->n) { //z
                gen[idx++] = this->z.at(i)[j];
            } else { //r
                if (this->r.at(i) == 1) {
                    gen[idx++] = true;
                } else {
                    gen[idx++] = false; //either 0 or 1 possible for phase
                }
            }
        }
        idx = 0;
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
    if (target > this->n || target < 0 ||
        control > this->n || control < 0) {
        return;
    }
    for (std::size_t i = 0U; i < this->n; ++i) {
        r[i] ^= (x[i][control] * z[i][target]) * (x[i][target] ^ z[i][control] ^ 1);
        x[i][target] ^= x[i][control];
        z[i][control] ^= z[i][target];
    }
}
void SatEncoder::QState::applyH(unsigned long target) {
    if (target > this->n || target < 0) {
        return;
    }
    for (std::size_t i = 0U; i < this->n; i++) {
        r[i] ^= x[i][target] * z[i][target];
        x[i][target] ^= z[i][target];
        z[i][target] = x[i][target] ^ z[i][target];
        x[i][target] ^= z[i][target];
    }
}
void SatEncoder::QState::applyS(unsigned long target) {
    if (target > this->n || target < 0) {
        return;
    }
    for (std::size_t i = 0U; i < this->n; ++i) {
        r[i] ^= x[i][target] * z[i][target];
        z[i][target] ^= x[i][target];
    }
}
void SatEncoder::QState::SetN(unsigned long n) {
    QState::n = n;
}
const std::vector<boost::dynamic_bitset<>>& SatEncoder::QState::GetX() const {
    return x;
}
void SatEncoder::QState::SetX(const std::vector<boost::dynamic_bitset<>>& x) {
    QState::x = x;
}
void SatEncoder::QState::SetZ(const std::vector<boost::dynamic_bitset<>>& z) {
    QState::z = z;
}
void SatEncoder::QState::SetR(const std::vector<int>& r) {
    QState::r = r;
}
const boost::uuids::uuid& SatEncoder::QState::GetPrevGenId() const {
    return prevGenId;
}
void SatEncoder::QState::SetPrevGenId(const boost::uuids::uuid& prev_gen_id) {
    prevGenId = prev_gen_id;
}
void SatEncoder::QState::printStateTableau() {
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
}
