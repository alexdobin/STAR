#ifndef __SAM_HEADER_H__
#define __SAM_HEADER_H__

#ifdef __cplusplus
extern "C" {
#endif

	void *sam_header_parse2(const char *headerText);
	void *sam_header_merge(int n, const void **dicts);
	void sam_header_free(void *header);
	char *sam_header_write(const void *headerDict);   // returns a newly allocated string

    /*
        // Usage example 
        const char *key, *val; 
        void *iter = sam_header_parse2(bam->header->text);
        while ( iter = sam_header_key_val(iter, "RG","ID","SM" &key,&val) ) printf("%s\t%s\n", key,val);
    */
    void *sam_header2key_val(void *iter, const char type[2], const char key_tag[2], const char value_tag[2], const char **key, const char **value);
	char **sam_header2list(const void *_dict, char type[2], char key_tag[2], int *_n);

    /*
        // Usage example
        int i, j, n;
        const char *tags[] = {"SN","LN","UR","M5",NULL}; 
        void *dict = sam_header_parse2(bam->header->text);
        char **tbl = sam_header2tbl_n(h->dict, "SQ", tags, &n);
        for (i=0; i<n; i++)
        {
            for (j=0; j<4; j++) 
                if ( tbl[4*i+j] ) printf("\t%s", tbl[4*i+j]); 
                else printf("-");
            printf("\n");
        }
        if (tbl) free(tbl);
     */
    char **sam_header2tbl_n(const void *dict, const char type[2], const char *tags[], int *n);

	void *sam_header2tbl(const void *dict, char type[2], char key_tag[2], char value_tag[2]);
	const char *sam_tbl_get(void *h, const char *key);
	int sam_tbl_size(void *h);
	void sam_tbl_destroy(void *h);

#ifdef __cplusplus
}
#endif

#endif
