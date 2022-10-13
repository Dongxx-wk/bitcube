//
//  main.cpp
//  bitcube
//
//  Created by Kosuke Shimizu on 2019/04/03.
//  Copyright © 2019 Kosuke Shimzu. All rights reserved.
//

#include "bitcube.hpp"

int main(int argc, const char * argv[]) {
    // insert code here...
    QDCT theQDCT(argv[1], argv[2], argv[3]);
    unsigned key[8] =\
    {0xc9e71c7c, 0xeee1a96c, 0x9b983b5f, 0x7442199e,\
        0xc261dfc4, 0x20141512, 0x3afc1a65, 0x617dbccd};
    theQDCT.input_key(key);

    theQDCT.bitcube_scrambling();

    theQDCT.write(argv[4], argv[5], argv[6]);
    return 0;
}
