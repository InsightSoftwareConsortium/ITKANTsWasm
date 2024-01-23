# the top-level README is used for describing this module, just
# re-used it for documentation here
get_filename_component(MY_CURRENT_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(READ "${MY_CURRENT_DIR}/README.rst" DOCUMENTATION)

# itk_module() defines the module dependencies in ANTsWasm
# ANTsWasm depends on ITKCommon
# The testing module in ANTsWasm depends on ITKTestKernel
# and ITKMetaIO(besides ANTsWasm and ITKCore)
# By convention those modules outside of ITK are not prefixed with
# ITK.

# define the dependencies of the include module and the tests
itk_module(ANTsWasm
  DEPENDS
    ITKCommon
    ITKTransform
    ITKRegistrationMethodsv4  # ImageRegistrationMethodv4 wraps DataObjectDecorator<CompositeTransform>
    ITKIOImageBase
    ITKImageLabel
    ITKBinaryMathematicalMorphology
    ITKTransformFactory
    ITKIOTransformBase
    ITKImageGrid
  TEST_DEPENDS
    ITKTestKernel
    ITKMetaIO
    ITKIONIFTI
    ITKIOTransformMatlab
    ITKIOTransformInsightLegacy
  DESCRIPTION
    "${DOCUMENTATION}"
  EXCLUDE_FROM_DEFAULT
  ENABLE_SHARED
)
