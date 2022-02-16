/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "CircuitOptimizer.hpp"
#include "algorithms/RandomCliffordCircuit.hpp"
#include "satencoder/SatEncoder.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <filesystem>
#include <locale>

class SatEncoderPTest: public testing::TestWithParam<std::string> {
protected:
    std::string            test_example_dir = "../../examples/";
    qc::QuantumComputation circuitOne{};
    qc::QuantumComputation circuitTwo{};

    void SetUp() override {
        circuitOne.import(test_example_dir + GetParam() + ".qasm");
        circuitTwo.import(test_example_dir + GetParam() + ".qasm");
    }
};
class SatEncoderTest: public testing::TestWithParam<std::string> {
protected:
    std::string            test_example_dir = "../../examples/";
    qc::QuantumComputation circuitOne{};
    qc::QuantumComputation circuitTwo{};
};
TEST_F(SatEncoderTest, CheckEqualWhenEqualRandomCircuits) {
    std::random_device        rd;
    std::mt19937              gen(rd());
    qc::RandomCliffordCircuit circOne(2, 1, gen());
    qc::CircuitOptimizer::flattenOperations(circOne);
    auto circTwo = circOne.clone();

    SatEncoder               sat_encoder;
    std::vector<std::string> inputs;
    std::string              filename{};
    bool                     result = sat_encoder.testEqual(circOne, circTwo, inputs, filename);
    EXPECT_EQ(result, true);
}

TEST_F(SatEncoderTest, CheckEqualWhenEqualRandomCircuitsWithInputs) {
    std::random_device        rd;
    std::mt19937              gen(rd());
    qc::RandomCliffordCircuit circOne(50, 10, gen());
    qc::CircuitOptimizer::flattenOperations(circOne);
    auto circTwo = circOne.clone();

    SatEncoder               sat_encoder;
    std::vector<std::string> inputs;
    inputs.emplace_back("ZX");
    inputs.emplace_back("ZZ");
    inputs.emplace_back("YZ");
    inputs.emplace_back("YY");
    inputs.emplace_back("XZ");
    std::string filename{};

    bool result = sat_encoder.testEqual(circOne, circTwo, inputs, filename);
    EXPECT_EQ(result, true);
}

INSTANTIATE_TEST_SUITE_P(SatEncoder, SatEncoderPTest,
                         testing::Values(
                                 "bell",
                                 "ghz",
                                 "simons"),
                         [](const testing::TestParamInfo<SatEncoderTest::ParamType>& info) {
                             std::string name = info.param;
                             std::replace(name.begin(), name.end(), '-', '_');
                             std::stringstream ss{};
                             ss << name;
                             return ss.str(); });
TEST_P(SatEncoderPTest, CheckEqualWhenEqualNoInputs) {
    SatEncoder               sat_encoder;
    std::vector<std::string> inputs;
    std::string              filename{};
    bool                     result = sat_encoder.testEqual(circuitOne, circuitTwo, inputs, filename);
    EXPECT_EQ(result, true);
}
TEST_P(SatEncoderPTest, CheckEqualWhenEqualTwoInputs) {
    SatEncoder               sat_encoder;
    std::vector<std::string> inputs;

    inputs.emplace_back("ZX");
    inputs.emplace_back("ZZ");
    std::string filename{};
    bool        result = sat_encoder.testEqual(circuitOne, circuitTwo, inputs, filename);
    EXPECT_EQ(result, true);
}

TEST_P(SatEncoderPTest, CheckEqualWhenNotEqualTwoInputs) {
    SatEncoder               sat_encoder;
    std::vector<std::string> inputs;

    inputs.emplace_back("ZX");
    inputs.emplace_back("ZZ");
    std::string filename{};
    circuitTwo.erase(circuitTwo.begin());
    bool result = sat_encoder.testEqual(circuitOne, circuitTwo, inputs, filename);
    EXPECT_EQ(result, true);
}

TEST_P(SatEncoderTest, SatWithoutInputs) {
    SatEncoder               sat_encoder;
    std::vector<std::string> inputs;
    std::string              file = "";
    sat_encoder.checkSatisfiability(circuitOne, inputs, file);
}
TEST_P(SatEncoderTest, SatWithInputs) {
    SatEncoder               sat_encoder;
    std::vector<std::string> inputs;
    inputs.emplace_back("XX");
    inputs.emplace_back("Zz");
    inputs.emplace_back("xZ");
    std::string file = "";
    sat_encoder.checkSatisfiability(circuitOne, inputs, file);
}

/* Benchmarking */
/*class SatEncoderBenchmarking: public testing::TestWithParam<std::string> {
protected:
    std::string test_example_dir = "../../examples/";
};
TEST_F(SatEncoderBenchmarking, GrowingNrOfQubitsForDepth10) {
    const size_t       depth         = 10;
    size_t             nrOfQubits    = 2;
    size_t             stepsize      = 10;
    size_t             maxNrOfQubits = 1000;
    std::random_device rd;
    std::ostringstream oss;

    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);
    oss << std::put_time(&tm, "%d-%m-%Y");
    auto          filename = oss.str();
    std::string   filep    = "/home/luca/Desktop/benchmark " + filename + ".json";
    std::ofstream outfile(filep, std::fstream::app);
    outfile << "{ \"benchmarks\" : [";
    outfile.close();
    for (; nrOfQubits < 30; nrOfQubits += stepsize) {
        SatEncoder               sat_encoder;
        std::vector<std::string> inputs;

        qc::RandomCliffordCircuit circOne(nrOfQubits, depth, rd());
        qc::CircuitOptimizer::flattenOperations(circOne);
        sat_encoder.checkSatisfiability(circOne, inputs, filep);
    }
    std::ofstream outfile2(filep, std::fstream::app);
    outfile2 << "]}";
    outfile2.close();
}

TEST_F(SatEncoderBenchmarking, GrowingNumberOfQubitsRandomFaults) {
    const size_t       depth         = 10;
    size_t             nrOfQubits    = 2;
    size_t             stepsize      = 10;
    size_t             maxNrOfQubits = 1000;
    std::random_device rd;
    std::ostringstream oss;
    std::mt19937       gen(rd());

    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);
    oss << std::put_time(&tm, "%d-%m-%Y");
    auto          timestamp = oss.str();
    std::string   filename    = "/home/luca/Desktop/benchmarkEC " + timestamp + ".json";
    std::ofstream outfile(filename, std::fstream::app);
    outfile << "{ \"benchmarks\" : [";
    outfile.close();
    for (; nrOfQubits < 30; nrOfQubits += stepsize) {
        SatEncoder               sat_encoder;
        std::vector<std::string> inputs;

        qc::RandomCliffordCircuit circOne(nrOfQubits, depth, gen());
        qc::CircuitOptimizer::flattenOperations(circOne);
        auto circTwo = circOne.clone();

        std::uniform_int_distribution<> distr(1, circTwo.size());
        int                             randomNumber;
        randomNumber = 1 + rand() % 2;
        if (randomNumber == 1) {
            circTwo.erase(circTwo.begin() + distr(gen));
        }
        sat_encoder.testEqual(circOne,circTwo, inputs, filename);
    }
    std::ofstream outfile2(filename, std::fstream::app);
    outfile2 << "]}";
    outfile2.close();

}*/