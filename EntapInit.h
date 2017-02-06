//
// Created by harta55 on 2/1/17.
//

#ifndef ENTAP_INITHANDLER_H
#define ENTAP_INITHANDLER_H

#include <string>


namespace entapInit {

    struct FtpFile {
        const char *filename;
        FILE *stream;
    };

    void download_file(std::string, std::string);

}

#endif //ENTAP_INITHANDLER_H
