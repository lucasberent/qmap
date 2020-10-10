//
// Created by Lukas Burgholzer on 28.02.20.
//

#include "Mapper.hpp"

void Mapper::initResults() {
	results.input_name = qc.getName();
	results.input_qubits = qc.getNqubits();
	results.architecture = architecture.getArchitectureName();
	results.calibration = architecture.getCalibrationName();
	results.output_name = qc.getName() + "_mapped";
	results.output_qubits = architecture.getNqubits();
	results.layeringStrategy = settings.layeringStrategy;
	results.initialLayoutStrategy = settings.initialLayoutStrategy;
	qcMapped.addQubitRegister(architecture.getNqubits());
}

Mapper::Mapper(const std::string& filename) {
	qc.import(filename);
	qc.stripIdleQubits(true);
	qubits.fill(DEFAULT_POSITION);
	locations.fill(DEFAULT_POSITION);
	fidelities.fill(INITIAL_FIDELITY);
}

/*
Mapper::Mapper(qc::QuantumComputation& quantumComputation) {
	qc = quantumComputation;
	qubits.fill(DEFAULT_POSITION);
	locations.fill(DEFAULT_POSITION);
}
 */

Mapper::Mapper(const std::string& filename, const std::string& cm_filename):
	Mapper(filename) {
	architecture = Architecture(cm_filename);
}

Mapper::Mapper(const std::string& filename, const std::string& cm_filename, const std::string& cal_filename):
	Mapper(filename, cm_filename) {
	architecture.loadCalibrationData(cal_filename);
}

Mapper::Mapper(const std::string& filename, unsigned short nQ, const CouplingMap& couplingMap):
	Mapper(filename) {
	architecture = Architecture(nQ, couplingMap);
}

Mapper::Mapper(const std::string& filename, unsigned short nQ, const CouplingMap& couplingMap,
		const std::vector<Architecture::CalibrationData>& calibrationData):
	Mapper(filename, nQ, couplingMap) {
	architecture.loadCalibrationData(calibrationData);
}

/*
Mapper::Mapper(qc::QuantumComputation& qc, const std::string& cm_filename):
	Mapper(qc) {
	architecture = Architecture(cm_filename);
}

Mapper::Mapper(qc::QuantumComputation& qc, const std::string& cm_filename, const std::string& cal_filename):
	Mapper(qc, cm_filename) {
	architecture.loadCalibrationData(cal_filename);
}

Mapper::Mapper(qc::QuantumComputation& qc, unsigned short nQ, const CouplingMap& couplingMap):
	Mapper(qc) {
	architecture = Architecture(nQ, couplingMap);
}

Mapper::Mapper(qc::QuantumComputation& qc, unsigned short nQ, const CouplingMap& couplingMap,
		const std::vector<Architecture::CalibrationData>& calibrationData):
	Mapper(qc, nQ, couplingMap) {
	architecture.loadCalibrationData(calibrationData);
}
 */

void Mapper::createLayers() {
	std::array<short, MAX_DEVICE_QUBITS> lastLayer{};
	lastLayer.fill(DEFAULT_POSITION);

	auto qubitsInLayer = std::set<unsigned short>{};

	for (auto& gate: qc) {
		if(!gate->isUnitary()) {
			throw QMAPException("Mapping not possible: circuit contains non-unitary operation: " + std::string(gate->getName()));
		}

		if (gate->getControls().size() > 1 || gate->getTargets().size() > 1) {
			throw QMAPException("Circuit contains gates with more than one control. Please make sure that the input circuit's gates are decomposed to the appropriate gate set!");
		}

		bool singleQubit = gate->getControls().empty();
		short control = -1;
		if (!singleQubit) {
			control = static_cast<short>(gate->getControls().at(0).qubit);
		}
		unsigned short target = gate->getTargets().at(0);
		size_t layer = 0;

		switch (settings.layeringStrategy) {
			case LayeringStrategy::IndividualGates:
			case LayeringStrategy::None:
				layers.emplace_back();
				layers.back().emplace_back(control, target, gate.get());
				break;
			case LayeringStrategy::DisjointQubits:
				if (singleQubit) {
					layer = lastLayer.at(target) + 1;
					lastLayer.at(target) = layer;
				} else {
					layer = std::max(lastLayer.at(control), lastLayer.at(target)) + 1;
					lastLayer.at(control) = lastLayer.at(target) = layer;
				}

				if (layers.size() <= layer) {
					layers.emplace_back();
				}
				layers.at(layer).emplace_back(control, target, gate.get());
				break;
			case LayeringStrategy::OddGates:
				if (results.input_gates % 2 == 0) {
					layers.emplace_back();
					layers.back().emplace_back(control, target, gate.get());
				} else {
					layers.back().emplace_back(control, target, gate.get());
				}
				break;
			case LayeringStrategy::QubitTriangle:
				if (layers.empty()) {
					layers.emplace_back();
				}

				if (singleQubit) {
					// single qubit gates can be added in any layer
					layers.back().emplace_back(control, target, gate.get());
				} else {
					qubitsInLayer.insert(control);
					qubitsInLayer.insert(target);

					if (qubitsInLayer.size() <= 3) {
						layers.back().emplace_back(control, target, gate.get());
					} else {
						layers.emplace_back();
						layers.back().emplace_back(control, target, gate.get());
						qubitsInLayer.clear();
						qubitsInLayer.insert(control);
						qubitsInLayer.insert(target);
					}
				}
				break;
		}

		if (singleQubit) {
			results.input_singlequbitgates++;
		} else {
			results.input_cnots++;
		}
		results.input_gates++;
	}
	results.input_layers = layers.size();
}

long Mapper::getNextLayer(size_t idx) {
	auto next = idx + 1;
	while (next < layers.size()) {
		for (const auto& gate: layers.at(next)) {
			if (!gate.singleQubit()) {
				return next;
			}
		}
		next++;
	}
	return -1;
}



