//
// Created by Lukas Burgholzer on 28.02.20.
//

#include "heuristic/HeuristicMapper.hpp"

void HeuristicMapper::map(const MappingSettings& ms) {
	settings = ms;
	if (settings.layeringStrategy == LayeringStrategy::OddGates || settings.layeringStrategy == LayeringStrategy::QubitTriangle) {
		std::cerr << "Layering strategy " << toString(settings.layeringStrategy) << " not suitable for heuristic mapper!" << std::endl;
		return;
	}
	auto start = std::chrono::high_resolution_clock::now();
	initResults();

	createLayers();
	if (ms.verbose) {
		printLayering(std::clog);
	}

	createInitialMapping();
	if (ms.verbose) {
		printLocations(std::clog );
		printQubits(std::clog );
	}

	unsigned long gateidx = 0;
	std::vector<unsigned long> gatesToAdjust{};
	for (size_t i = 0; i < layers.size(); ++i) {
		Node result = AstarMap(i);

		qubits = result.qubits;
		locations = result.locations;

		// initial layer needs no swaps
		if(i != 0) {
			for (const auto& swaps: result.swaps) {
				for (const auto& swap: swaps) {
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), std::vector<qc::Control>{ }, swap.first, swap.second, qc::SWAP);
					std::swap(qcMapped.outputPermutation.at(swap.first), qcMapped.outputPermutation.at(swap.second));
					results.output_swaps++;
					if (architecture.bidirectional()) {
						results.output_gates += GATES_OF_BIDIRECTIONAL_SWAP;
					} else {
						results.output_gates += GATES_OF_UNIDIRECTIONAL_SWAP;
					}
					gateidx++;
				}
			}
		}

		// add gates of the layer to circuit
		for(const auto& gate: layers.at(i)) {
			auto op = dynamic_cast<qc::StandardOperation*>(gate.op);
			if (!op) {
				throw QMAPException("Cast to StandardOperation not possible during mapping. Check that circuit contains only StandardOperations");
			}

			if (gate.singleQubit()) {
				if (locations.at(gate.target) == DEFAULT_POSITION) {
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(),
					                                             gate.target,
					                                             op->getType(),
					                                             op->getParameter().at(0),
					                                             op->getParameter().at(1),
					                                             op->getParameter().at(2));
					gatesToAdjust.push_back(gateidx);
					results.output_gates++;
					results.output_singlequbitgates++;
					gateidx++;
				} else {
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(),
							locations.at(gate.target),
							op->getType(),
							op->getParameter().at(0),
							op->getParameter().at(1),
							op->getParameter().at(2));
					results.output_gates++;
					results.output_singlequbitgates++;
					gateidx++;
				}
			} else {
				Edge cnot = {locations.at(gate.control), locations.at(gate.target)};
				if (architecture.getCouplingMap().find(cnot) == architecture.getCouplingMap().end()) {
					Edge reverse = {cnot.second, cnot.first};
					if (architecture.getCouplingMap().find(reverse) == architecture.getCouplingMap().end()) {
						throw QMAPException("Invalid CNOT: " + std::to_string(reverse.first) + "-" + std::to_string(reverse.second));
					}
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), reverse.first, qc::H);
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), reverse.second, qc::H);
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), qc::Control(reverse.first), reverse.second, qc::X);
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), reverse.second, qc::H);
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(), reverse.first, qc::H);

					results.output_direction_reverse++;
					results.output_cnots++;
					results.output_gates += 5;
					gateidx += 5;
				} else {
					qcMapped.emplace_back<qc::StandardOperation>(qcMapped.getNqubits(),
					                                             qc::Control(cnot.first), cnot.second, qc::X);
					results.output_cnots++;
					results.output_gates++;
					gateidx++;
				}
			}

		}
	}

	// fix single qubit gates
	if (!gatesToAdjust.empty()) {
		gateidx--; // index of last operation
		for(auto it = qcMapped.rbegin(); it != qcMapped.rend(); ++it, --gateidx) {
			auto op = dynamic_cast<qc::StandardOperation*>(it->get());
			if (!op) {
				throw QMAPException("Cast to StandardOperation not possible during mapping. Check that circuit contains only StandardOperations");
			}
			if (op->getType() == qc::SWAP) {
				short q0 = qubits.at(op->getTargets().at(0));
				short q1 = qubits.at(op->getTargets().at(1));
				qubits.at(op->getTargets().at(0)) = q1;
				qubits.at(op->getTargets().at(1)) = q0;

				if (q0 != DEFAULT_POSITION) {
					locations.at(q0) = op->getTargets().at(1);
				}
				if (q1 != DEFAULT_POSITION) {
					locations.at(q1) = op->getTargets().at(0);
				}
			}
			if (gatesToAdjust.back() == gateidx) {
				gatesToAdjust.pop_back();
				auto target = op->getTargets().at(0);
				auto targetLocation = locations.at(target);

				if (targetLocation == -1) {
					// qubit only occurs in single qubit gates, can be mapped to an arbitrary free qubit
					unsigned short loc = 0;
					while(qubits.at(loc) != DEFAULT_POSITION) {
						++loc;
					}
					locations.at(target) = loc;
					op->setTargets({loc});
					qcMapped.initialLayout.at(target) = loc;
					qc.outputPermutation.at(target) = loc;
				} else {
					op->setTargets({ static_cast<unsigned short>(targetLocation)});
				}
			}
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = end - start;
	results.time = diff.count();
	results.timeout = false;
}

void HeuristicMapper::initResults() {
	Mapper::initResults();
	results.method = Method::Heuristic;
}

void HeuristicMapper::createInitialMapping() {
	if (layers.empty())
		return;

	switch (settings.initialLayoutStrategy) {
		case InitialLayoutStrategy::Identity:
			for (unsigned short i = 0; i < qc.getNqubits(); ++i) {
				locations.at(i) = i;
				qubits.at(i) = i;
			}
			break;
		case InitialLayoutStrategy::Static:
			for (const auto& gate: layers.at(0)) {
				if (gate.singleQubit())
					continue;

				for (const auto& edge: architecture.getCouplingMap()) {
					if (qubits.at(edge.first) == DEFAULT_POSITION && qubits.at(edge.second) == DEFAULT_POSITION) {
						qubits.at(edge.first) = gate.control;
						qubits.at(edge.second) = gate.target;
						locations.at(gate.control) = edge.first;
						locations.at(gate.target) = edge.second;
					}
				}
			}

			// assign remaining logical qubits
			for (unsigned short i = 0; i < qc.getNqubits(); ++i) {
				if (locations.at(i) == DEFAULT_POSITION) {
					for (unsigned short j=0; j<architecture.getNqubits(); ++j) {
						if (qubits.at(j) == DEFAULT_POSITION) {
							locations.at(i) = j;
							qubits.at(j) = i;
							break;
						}
					}
				}
			}

			break;
		case InitialLayoutStrategy::Dynamic:
		case InitialLayoutStrategy::None:
			// nothing to be done here
			break;

		// TODO: Design strategy that maps most used qubit to most connected qubits on architecture
	}
}

void HeuristicMapper::mapUnmappedGates(long layer, HeuristicMapper::Node& node, std::vector<unsigned short>& consideredQubits) {
	for (const auto& gate: layers.at(layer)) {
		if (gate.singleQubit())
			continue;

		consideredQubits.push_back(gate.control);
		consideredQubits.push_back(gate.target);

		auto controlLocation = locations.at(gate.control);
		auto targetLocation = locations.at(gate.target);

		if (controlLocation == DEFAULT_POSITION && targetLocation == DEFAULT_POSITION) {
			std::set<Edge> possibleEdges{};
			for (const auto& edge: architecture.getCouplingMap()) {
				if (qubits.at(edge.first) == DEFAULT_POSITION && qubits.at(edge.second) == DEFAULT_POSITION) {
					possibleEdges.emplace(edge);
				}
			}

			if (possibleEdges.empty()) {
				throw QMAPException("Could not map logical qubits to physical qubits. No suitable edge found.");
			}

			// TODO: Consider fidelity here if available. The best available edge should be chosen
			auto& chosenEdge = *possibleEdges.begin();
			locations.at(gate.control) = chosenEdge.first;
			locations.at(gate.target) = chosenEdge.second;
			qubits.at(chosenEdge.first) = gate.control;
			qubits.at(chosenEdge.second) = gate.target;
			qcMapped.initialLayout.at(chosenEdge.first) = gate.control;
			qcMapped.initialLayout.at(chosenEdge.second) = gate.target;
			qcMapped.outputPermutation.at(chosenEdge.first) = gate.control;
			qcMapped.outputPermutation.at(chosenEdge.second) = gate.target;
		} else if (controlLocation == DEFAULT_POSITION) {
			mapToMinDistance(gate.target, gate.control);
		} else if (targetLocation == DEFAULT_POSITION) {
			mapToMinDistance(gate.control, gate.target);
		}
		node.costHeur = std::max(node.costHeur, distanceOnArchitectureOfLogicalQubits(gate.control, gate.target));
	}

	node.done = node.costHeur <= COST_DIRECTION_REVERSE;
}

void HeuristicMapper::mapToMinDistance(unsigned short source, unsigned short target) {
	auto min = std::numeric_limits<double>::max();
	unsigned short pos = DEFAULT_POSITION;
	for (int i = 0; i < architecture.getNqubits(); ++i) {
		if (qubits.at(i) == DEFAULT_POSITION) {
			// TODO: Consider fidelity here if available
			auto distance = distanceOnArchitectureOfPhysicalQubits(locations.at(source), i);
			if (distance < min) {
				min = distance;
				pos = i;
			}
		}
	}
	qubits.at(pos) = target;
	locations.at(target) = pos;
	qcMapped.initialLayout.at(pos) = target;
	qcMapped.outputPermutation.at(pos) = target;
}

HeuristicMapper::Node HeuristicMapper::AstarMap(long layer) {
	std::vector<unsigned short> consideredQubits{};
	Node node{};

	mapUnmappedGates(layer, node, consideredQubits);

	node.locations = locations;
	node.qubits = qubits;

	nodes.push(node);
	std::vector<Edge> swaps{consideredQubits.size()};

	while (!nodes.top().done) {
		Node current = nodes.top();
		nodes.pop();

		expandNode(consideredQubits, 0, swaps, 0, current, layer);
	}

	Node result = nodes.top();
	nodes.pop();

	while (!nodes.empty()) {
		nodes.pop();
	}

	return result;

}

// TODO: replace with one swap per expand method
void HeuristicMapper::expandNode(const std::vector<unsigned short>& consideredQubits, unsigned short qubit, std::vector<Edge>& swaps, unsigned short nswaps, HeuristicMapper::Node& node, long layer) {
	auto& currentLayer = layers.at(layer);

	if (qubit == consideredQubits.size()) {
		if (nswaps == 0)
			return;

		Node n = Node(node.qubits, node.locations, node.swaps);
		n.nswaps = node.nswaps + nswaps;
		if (architecture.bidirectional()) {
			n.costFixed = node.costFixed + COST_BIDIRECTIONAL_SWAP * swaps.size();
		} else {
			n.costFixed = node.costFixed + COST_UNIDIRECTIONAL_SWAP * swaps.size();
		}

		n.swaps.emplace_back();
		for (int i = 0; i < nswaps; ++i) {
			n.applySWAP(swaps.at(i));
		}
		n.done = true;

		for (const auto& gate: currentLayer) {
			if (gate.singleQubit())
				continue;

			n.updateHeuristicCost(architecture, gate, settings.admissibleHeuristic);
			n.checkUnfinished(architecture, gate);
		}

		// calculate heuristics for the cost of the following layers
		if (settings.lookahead) {
			lookahead(getNextLayer(layer), n);
		}

		nodes.push(n);
	} else {
		expandNode(consideredQubits, qubit+1, swaps, nswaps, node, layer);
		for(const auto& edge: architecture.getCouplingMap()) {
			if (edge.first == node.locations.at(consideredQubits.at(qubit)) || edge.second == node.locations.at(consideredQubits.at(qubit))) {
				if (usedDeviceQubits.count(edge.first) == 0 && usedDeviceQubits.count(edge.second) == 0) {
					usedDeviceQubits.insert(edge.first);
					usedDeviceQubits.insert(edge.second);
					swaps.at(nswaps) = edge;
					expandNode(consideredQubits, qubit+1, swaps, nswaps+1, node, layer);
					usedDeviceQubits.erase(edge.first);
					usedDeviceQubits.erase(edge.second);
				}
			}
		}
	}
}

void HeuristicMapper::lookahead(long layer, HeuristicMapper::Node& node) {
	auto nextLayer = layer;
	double factor = settings.firstLookaheadFactor;
	for (int i = 0; i < settings.nrLookaheads; ++i) {
		if (nextLayer == -1) {
			break;
		}

		double penalty = 0.;
		for (const auto& gate: layers.at(nextLayer)) {
			if (gate.singleQubit())
				continue;

			auto loc1 = node.locations.at(gate.control);
			auto loc2 = node.locations.at(gate.target);
			if (loc1 == DEFAULT_POSITION && loc2 == DEFAULT_POSITION) {
				// no penalty
			} else if (loc1 == DEFAULT_POSITION) {
				auto min = std::numeric_limits<double>::max();
				for (int j = 0; j < architecture.getNqubits(); ++j) {
					if (node.qubits.at(j) == DEFAULT_POSITION) {
						// TODO: Consider fidelity here if available
						min = std::min(min, distanceOnArchitectureOfPhysicalQubits(j, node.locations.at(gate.target)));
					}
				}
				penalty = heuristicCost(penalty, min);
			} else if (loc2 == DEFAULT_POSITION) {
				auto min = std::numeric_limits<double>::max();
				for (int j = 0; j < architecture.getNqubits(); ++j) {
					if (node.qubits.at(j) == DEFAULT_POSITION) {
						// TODO: Consider fidelity here if available
						min = std::min(min, distanceOnArchitectureOfPhysicalQubits(node.locations.at(gate.target), j));
					}
				}
				penalty = heuristicCost(penalty, min);
			} else {
				auto cost = architecture.distance(node.locations.at(gate.control), node.locations.at(gate.target));
				penalty = heuristicCost(penalty, cost);
			}
		}

		node.lookaheadPenalty += factor * penalty;
		factor *= settings.lookaheadFactor;
		nextLayer = getNextLayer(nextLayer); // TODO: consider single qubits here for better fidelity lookahead
	}

}