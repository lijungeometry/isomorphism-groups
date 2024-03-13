/* lib.d file reduce.h */
extern void diff_reduction PARMS((twoafsa* DIFF,word * wp));
extern int *** reduction_rules PARMS((word ** rulewords,int nrules));
extern void wa_reduction PARMS((afsa* WA,int *** rules,word * wp));
extern void clear_reduction_rules PARMS((int *** rules));


