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

#ifndef __UNITTESTSUMMARY_H
#define __UNITTESTSUMMARY_H

#include "Output.hpp"
#include "UnitTestResult.hpp"

namespace MainApplication
{
    class UnitTestSummary;

    class UnitTestSummary
    {
        protected:
            static map<unsigned int, UnitTestResult> testPerformed; ///< unitTest results

            static unsigned int numberTotalTests;
            static unsigned int numberSuccessTests;
            static unsigned int numberErrorTests;
        public:
            static const unsigned int& NumberTotalTests() { return numberTotalTests; }
            static const unsigned int& NumberSuccessTests() { return numberSuccessTests; }
            static const unsigned int& NumberErrorTests() { return numberErrorTests; }

            /// Check the UnitTest result and create a new UnitTestResult
            /// @return the UnitTest UId
            static unsigned int CheckUnitTest(const bool& logicResult, const string& message);
            /// Check the UnitTest result and create a new UnitTestResult
            /// @return the UnitTest UId
            static unsigned int CheckUnitTest(const Output::ExitCodes& logicResult, const string& message);

            /// Reset the Summary
            static void Reset();
            /// Print the UnitTest summary
            static void PrintSummary();
            /// Export the UnitTest results to filePath
            static void ExportTestResult(const string& filePath, const char separator = ';');
    };
}

#endif // __UNITTESTSUMMARY_H
