/* One-time helper: generates binary blobs for each of the 50 metagenome
   training models. Each blob is a raw dump of struct _training. */

#include <stdio.h>
#include <string.h>
#include "training.h"

typedef void (*init_fn)(struct _training *);

int main() {
    struct _training tinf;
    char path[256];
    FILE *fp;
    int i;

    init_fn fns[50] = {
        initialize_metagenome_0,  initialize_metagenome_1,
        initialize_metagenome_2,  initialize_metagenome_3,
        initialize_metagenome_4,  initialize_metagenome_5,
        initialize_metagenome_6,  initialize_metagenome_7,
        initialize_metagenome_8,  initialize_metagenome_9,
        initialize_metagenome_10, initialize_metagenome_11,
        initialize_metagenome_12, initialize_metagenome_13,
        initialize_metagenome_14, initialize_metagenome_15,
        initialize_metagenome_16, initialize_metagenome_17,
        initialize_metagenome_18, initialize_metagenome_19,
        initialize_metagenome_20, initialize_metagenome_21,
        initialize_metagenome_22, initialize_metagenome_23,
        initialize_metagenome_24, initialize_metagenome_25,
        initialize_metagenome_26, initialize_metagenome_27,
        initialize_metagenome_28, initialize_metagenome_29,
        initialize_metagenome_30, initialize_metagenome_31,
        initialize_metagenome_32, initialize_metagenome_33,
        initialize_metagenome_34, initialize_metagenome_35,
        initialize_metagenome_36, initialize_metagenome_37,
        initialize_metagenome_38, initialize_metagenome_39,
        initialize_metagenome_40, initialize_metagenome_41,
        initialize_metagenome_42, initialize_metagenome_43,
        initialize_metagenome_44, initialize_metagenome_45,
        initialize_metagenome_46, initialize_metagenome_47,
        initialize_metagenome_48, initialize_metagenome_49,
    };

    printf("sizeof(struct _training) = %zu\n", sizeof(struct _training));

    for (i = 0; i < 50; i++) {
        memset(&tinf, 0, sizeof(tinf));
        fns[i](&tinf);
        sprintf(path, "data/metagenome_%d.bin", i);
        fp = fopen(path, "wb");
        if (!fp) {
            fprintf(stderr, "Error: cannot open %s\n", path);
            return 1;
        }
        fwrite(&tinf, sizeof(struct _training), 1, fp);
        fclose(fp);
        printf("Wrote %s (%zu bytes)\n", path, sizeof(struct _training));
    }
    return 0;
}
