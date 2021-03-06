MACRO (CHANGE_FILE_EXTENSION FILE_EXT1 FILE_EXT2 OUTVAR LIST)
   SET(BLA)
   IF (${FILE_EXT1} MATCHES "^[*][.]+.*$")
     STRING(REGEX REPLACE "^[*]+([.].*)$" "\\1" FILE_EXT1_NEW ${FILE_EXT1}) 
   ENDIF  (${FILE_EXT1} MATCHES "^[*][.]+.*$")
   IF (${FILE_EXT2} MATCHES "^[*][.]+.*$")
     STRING(REGEX REPLACE "^[*]+([.].*)" "\\1" FILE_EXT2_NEW ${FILE_EXT2}) 
   ENDIF  (${FILE_EXT2} MATCHES "^[*][.]+.*$")
   foreach (_current_FILE ${LIST})
     STRING(REGEX REPLACE "^(.*)${FILE_EXT1_NEW}$" "\\1${FILE_EXT2_NEW}" test ${_current_FILE})
     SET (BLA ${BLA} ${test})
   endforeach (_current_FILE ${ARGN})
   SET (${OUTVAR} ${BLA})
ENDMACRO (CHANGE_FILE_EXTENSION)