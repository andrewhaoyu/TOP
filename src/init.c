#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern void CompleteCasesScoreTestSupport(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void convert(void *, void *, void *);
extern void EMStep(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void EMStepScoreSupport(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void EMStepScoreTestSupportMixedModel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Mvpoly_complete(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void OneStepMLE(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ScoreTest(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ScoreTestMixedModel(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"CompleteCasesScoreTestSupport",    (DL_FUNC) &CompleteCasesScoreTestSupport,    22},
    {"convert",                          (DL_FUNC) &convert,                           3},
    {"EMStep",                           (DL_FUNC) &EMStep,                           22},
    {"EMStepScoreSupport",               (DL_FUNC) &EMStepScoreSupport,               26},
    {"EMStepScoreTestSupportMixedModel", (DL_FUNC) &EMStepScoreTestSupportMixedModel, 24},
    {"Mvpoly_complete",                  (DL_FUNC) &Mvpoly_complete,                  18},
    {"OneStepMLE",                       (DL_FUNC) &OneStepMLE,                       21},
    {"ScoreTest",                        (DL_FUNC) &ScoreTest,                        20},
    {"ScoreTestMixedModel",              (DL_FUNC) &ScoreTestMixedModel,              19},
    {NULL, NULL, 0}
};

void attribute_visible R_init_TOP(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
