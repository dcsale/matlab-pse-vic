/*
 * stringify.c
 *   Turns input from stdin into a C string definition written on
 *   stdout.  Adds quotes, newline characters, and backslashes as
 *   needed.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */
#include <stdio.h>
#include <stdlib.h>


void stringify(const char* name)
{
    char line[512];
    char* p;

    printf("/*\n"
           " * Auto-generated by stringify\n"
           " */\n\n");
    printf("const char* %s =", name);
    while (fgets(line, sizeof(line), stdin)) {
        printf("\n  \"");
        for (p = line; *p; ++p) {
            if (*p == '"' || *p == '\\')
                putc('\\', stdout);
            if (*p != '\n' && *p != '\r')
                putc(*p, stdout);
        }
        printf("\\n\"");
    }
    printf(";\n\n");
}


int main(int argc, char** argv)
{
    if (argc != 2) {
        fprintf(stderr, "String name required\n");
        exit(-1);
    }
    stringify(argv[1]);
    return 0;
}
