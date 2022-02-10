/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "algorithms/RandomCliffordCircuit.hpp"
#include "satencoder/SatEncoder.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <locale>

class SatEncoderTest: public testing::TestWithParam<std::string> {
protected:
    std::string test_example_dir      = "./examples/";
    qc::QuantumComputation circuitOne{};
    qc::QuantumComputation circuitTwo{};

    void SetUp() override {
    }
};
TEST_F(SatEncoderTest, CheckEqualWhenEqual) {

    circuitOne.import(test_example_dir + "bell.qasm");
    circuitTwo.import(test_example_dir + "bell.qasm");
    SatEncoder               sat_encoder;
    std::vector<std::string> inputs;

    bool result = sat_encoder.testEqual(circuitOne, circuitTwo, inputs);
    std::cout << "Equal: " << result;
    EXPECT_EQ(result, true);
}

TEST_F(SatEncoderTest, SatTwoQubitBellCircuit) {
    circuitOne.import(test_example_dir + "bell.qasm");
    SatEncoder               sat_encoder;
    std::vector<std::string> inputs;

    sat_encoder.checkSatisfiability(circuitOne, inputs);
}
/*

TEST_P(SatEncoderTest, SatRandomCliffordCircuitGrowingSize) {
    SatEncoder               enc;
    std::vector<std::string> input; //empty == |0...0>
    size_t                   maxNrOfQubits = 100;
    long                     stepSize      = 5;
    long                     depth         = 10;

    for (size_t i = 0; i < maxNrOfQubits; i += stepSize) {
        qc::RandomCliffordCircuit randomCircuit(i, 10, 0);
        enc.checkSatisfiability(randomCircuit, input);
        randomCircuit.printStatistics(std::cout);
    }
}

*/
