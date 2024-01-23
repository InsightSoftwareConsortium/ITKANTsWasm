itk_module_test()

set(ANTsWasmTests
  itkANTSRegistrationTest.cxx
  itkANTSRegistrationBasicTests.cxx
  )

CreateTestDriver(ANTsWasm "${ANTsWasm-Test_LIBRARIES}" "${ANTsWasmTests}")


itk_add_test(NAME itkANTSRegistrationBasicTests
  COMMAND ANTsWasmTestDriver
  itkANTSRegistrationBasicTests  # parameterless test
  )

itk_add_test(NAME itkANTSRegistrationTest
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_AffineRotationNoMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineRotationNoMasks.result.nii.gz
    --compareIntensityTolerance 7
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 777
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/rotation.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineRotationNoMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineRotationNoMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_AffineRotationNoMasks.mat}  # initial transform
  )
