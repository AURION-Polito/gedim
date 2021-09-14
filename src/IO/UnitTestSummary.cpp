// Copyright (C) 2020 Vicini Fabio
//
// This file is part of the dissertation of the author.
//
// This is a free program: you can redistribute it and/or modify
// it under the terms of the author.
//
// This program is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.
//
// Modified by Vicini Fabio 2020
//
// First added:  2020-01-21

#include "UnitTestSummary.hpp"
#include "MpiParallelEnvironment.hpp"

using namespace GeDiM;

namespace MainApplication
{
    // ***************************************************************************
    map<unsigned int, UnitTestResult> UnitTestSummary::testPerformed = map<unsigned int, UnitTestResult>();
    unsigned int UnitTestSummary::numberTotalTests = 0;
    unsigned int UnitTestSummary::numberSuccessTests = 0;
    unsigned int UnitTestSummary::numberErrorTests = 0;
    // ***************************************************************************
    unsigned int UnitTestSummary::CheckUnitTest(const bool& logicResult, const string& message)
    {
        if (MpiParallelEnvironment::Process().Rank() > 0)
            return 0;

        numberTotalTests++;

        auto unitTestResult = UnitTestResult();
        unitTestResult.UId = numberTotalTests;
        unitTestResult.Message = message;
        unitTestResult.Result = logicResult ? Output::Success : Output::GenericError;

        if (logicResult)
            numberSuccessTests++;
        else
            numberErrorTests++;

        testPerformed.insert(pair<unsigned int, UnitTestResult>(unitTestResult.UId, unitTestResult));

        return unitTestResult.UId;
    }
    // ***************************************************************************
    unsigned int UnitTestSummary::CheckUnitTest(const Output::ExitCodes& logicResult, const string &message)
    {
        if (MpiParallelEnvironment::Process().Rank() > 0)
            return 0;

        numberTotalTests++;

        auto unitTestResult = UnitTestResult();
        unitTestResult.UId = numberTotalTests;
        unitTestResult.Message = message;
        unitTestResult.Result = logicResult;

        switch (logicResult)
        {
            case Output::Success:
                numberSuccessTests++;
                break;
            default:
                numberErrorTests++;
                break;
        }

        testPerformed.insert(pair<unsigned int, UnitTestResult>(unitTestResult.UId, unitTestResult));

        return unitTestResult.UId;
    }
    // ***************************************************************************
    void UnitTestSummary::Reset()
    {
        if (MpiParallelEnvironment::Process().Rank() > 0)
            return;

        numberTotalTests = 0;
        numberSuccessTests = 0;
        numberErrorTests = 0;
        testPerformed.clear();
    }
    // ***************************************************************************
    void UnitTestSummary::PrintSummary()
    {
        if (MpiParallelEnvironment::Process().Rank() > 0)
            return;

        Output::PrintGenericMessage("Total Tests done: %d, Number Success Tests: %d, Number Error Tests: %d", true, numberTotalTests, numberSuccessTests, numberErrorTests);
    }
    // ***************************************************************************
    void UnitTestSummary::ExportTestResult(const string& filePath, const char separator)
    {
        if (MpiParallelEnvironment::Process().Rank() > 0)
            return;

        ofstream outFile(filePath.c_str());

        outFile<< "UId"<< separator;
        outFile<< "Result"<< separator;
        outFile<< "Description"<< separator;
        outFile<< "Total"<< separator;
        outFile<< "Successfull"<< separator;
        outFile<< "Error"<< endl;

        for (map<unsigned int, UnitTestResult>::iterator t = testPerformed.begin(); t != testPerformed.end(); t++ )
        {
            const UnitTestResult& testResult = t->second;

            outFile << testResult.UId<< separator;
            outFile << testResult.Result<< separator;
            outFile << Output::ExitCodeToString(testResult.Result, false)<< separator;
            outFile << numberTotalTests<< separator;
            outFile << numberSuccessTests<< separator;
            outFile << numberErrorTests<< endl;
        }
    }
	// ***************************************************************************
}
