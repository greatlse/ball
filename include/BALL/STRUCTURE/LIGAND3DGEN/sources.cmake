#########################################
# by pbrach                             #
#########################################

### list all filenames of the directory here ###
SET(GROUP LIGAND3DGEN)
 
SET(DIRECTORY include/BALL/STRUCTURE/LIG3DGEN)

FILE(GLOB HEADERS_LIST "${DIRECTORY}/*.h")

### add full path to the filenames ###
SET(HEADERS)
FOREACH(i ${HEADERS_LIST})
    ### make sure we do not have absolute paths flying around...
    GET_FILENAME_COMPONENT(i ${i} NAME)
    LIST(APPEND HEADERS ${DIRECTORY}/${i})
ENDFOREACH()

### pass source file list to the upper instance ###
SET(BALL_headers ${BALL_headers} ${HEADERS})

### source group definition ###
STRING(REGEX REPLACE "/" "\\\\" S_GROUP ${GROUP})
SOURCE_GROUP("Header Files\\\\${S_GROUP}" FILES ${HEADERS})