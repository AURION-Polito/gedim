#!/bin/bash

FILE=${1}
LICENCE_TEST="_LICENSE_HEADER_"
LICENSE_STRING="\/\/ _LICENSE_HEADER_\n\/\/\n\/\/ Copyright (C) 2019 - 2025.\n\/\/ Terms register on the GPL-3.0 license.\n\/\/\n\/\/ This file can be redistributed and/or modified under the license terms.\n\/\/\n\/\/ See top level LICENSE file for more details.\n\/\/\n\/\/ This file can be used citing references in CITATION.cff file.\n"

grep -qF "${LICENCE_TEST}" ${FILE} || sed -i "1i ${LICENSE_STRING}" ${FILE}
sed -i "1,11d" ${FILE}
grep -qF "${LICENCE_TEST}" ${FILE} || sed -i "1i ${LICENSE_STRING}" ${FILE}
