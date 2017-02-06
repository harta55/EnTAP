//
// Created by harta55 on 2/1/17.
//

#include <cstdio>
#include "EntapInit.h"
#include <curl.h>

namespace entapInit {

    size_t my_fwrite(void *buffer, size_t size, size_t nmemb, void *stream) {
        struct FtpFile *out = (struct FtpFile *) stream;
        if (out && !out -> stream) {
            // open file to write
            out -> stream = fopen(out->filename, "wb");
            if (!out->stream) return 1; // can't open file to write
        }
        return fwrite(buffer, size, nmemb, out->stream);
    }

    void download_file(std::string, std::string) {
        CURL *curl;
        CURLcode res;
        struct FtpFile ftpfile={
                "yourfile.bin", /* name to store the file as if successful */
                NULL
        };

        curl_global_init(CURL_GLOBAL_DEFAULT);

        curl = curl_easy_init();
        if(curl) {
            /*
             * You better replace the URL with one that works! Note that we use an
             * FTP:// URL with standard explicit FTPS. You can also do FTPS:// URLs if
             * you want to do the rarer kind of transfers: implicit.
             */
            curl_easy_setopt(curl, CURLOPT_URL,
                             "ftp://user@server/home/user/file.txt");
            /* Define our callback to get called when there's data to be written */
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, my_fwrite);
            /* Set a pointer to our struct to pass to the callback */
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, &ftpfile);

            /* We activate SSL and we require it for both control and data */
            curl_easy_setopt(curl, CURLOPT_USE_SSL, CURLUSESSL_ALL);

            /* Switch on full protocol/debug output */
            curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);

            res = curl_easy_perform(curl);

            /* always cleanup */
            curl_easy_cleanup(curl);

            if(CURLE_OK != res) {
                /* we failed */
                fprintf(stderr, "curl told us %d\n", res);
            }
        }

        if(ftpfile.stream)
            fclose(ftpfile.stream); /* close the local file */

        curl_global_cleanup();
    }
}
