// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME kkpDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "kkp.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace {
  void TriggerDictionaryInitialization_kkpDict_Impl() {
    static const char* headers[] = {
"kkp.h",
nullptr
    };
    static const char* includePaths[] = {
"/Users/gijsvanweelden/alice/sw/osx_x86-64/ROOT/v6-26-04-patches-alice1-local1/include/",
"/Users/gijsvanweelden/Documents/Fragmentation/brick/kkp/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "kkpDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "kkpDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "kkp.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"kkp", payloadCode, "@",
"kkp_func", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("kkpDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_kkpDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_kkpDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_kkpDict() {
  TriggerDictionaryInitialization_kkpDict_Impl();
}
