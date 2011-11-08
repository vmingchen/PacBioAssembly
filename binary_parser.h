/*
 * ===========================================================================
 *
 *       Filename:  binary_parser.h
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/08/2011 11:59:58 AM
 *
 *    Description:  
 *
 *       Revision:  none
 *
 * ===========================================================================
 */
#ifdef BINARY_PARSER_H
#define BINARY_PARSER_H

class binary_parser {
public:
//    // create binary file from text file 
//    void create(const char *in_path, const char *out_path);
    // read binary file
    bool parse(const char *path);
    // is empty or not
    bool empty();
    // read a record from the binary
    const unsigned* read(unsigned *plen) const;
    // convert text file to binary file, return length of result
    unsigned text2bin(const char *ptext, unsigned *pbin, unsigned blen);
    unsigned bin2text(const unsigned *pbin, char *ptext, unsigned tlen);
};

#endif
