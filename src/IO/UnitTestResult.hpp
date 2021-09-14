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

#ifndef __UNITTESTRESULT_H
#define __UNITTESTRESULT_H

#include "Output.hpp"

namespace MainApplication
{
    struct UnitTestResult;

    struct UnitTestResult
    {
        unsigned int UId; ///< Unique Id of the UnitTest
        Output::ExitCodes Result; ///< Result of the UnitTest
        string Message; ///< Message of the UnitTest
    };
}

#endif // __UNITTESTRESULT_H
