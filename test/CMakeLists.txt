itk_module_test()

set(ANTsWasmTests
  itkANTSRegistrationTest.cxx
  itkANTSRegistrationBasicTests.cxx
  )

CreateTestDriver(ANTsWasm "${ANTsWasm-Test_LIBRARIES}" "${ANTsWasmTests}")


itk_add_test(NAME itkANTSRegistrationBasicTests
  COMMAND ANTsWasmTestDriver
  itkANTSRegistrationBasicTests ${ITK_TEST_OUTPUT_DIR}
  )

itk_add_test(NAME antsRegistrationTest_AffineScaleMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_AffineScaleMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineScaleMasks.result.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/scale.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineScaleMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineScaleMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_AffineScaleMasks.mat}  # initial transform
    DATA{Input/test_mask.nii.gz}  # fixedMask
    DATA{Input/scale.test_mask.nii.gz}  # movingMask
    0.001  # (minimum?) StepLength
    Affine
    Mattes
    0.10  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    1x1x1  # shrinkFactors
    0.1x0.1x0.1 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_AffineScaleNoMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_AffineScaleNoMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineScaleNoMasks.result.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/scale.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineScaleNoMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineScaleNoMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_AffineScaleMasks.mat}  # initial transform
    none  # fixedMask
    none  # movingMask
    0.001  # (minimum?) StepLength
    Affine
    Mattes
    0.10  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    1x1x1  # shrinkFactors
    0.1x0.1x0.1 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_AffineRotationMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_AffineRotationMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineRotationMasks.result.nii.gz
    --compareIntensityTolerance 25
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/rotation.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineRotationMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineRotationMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_AffineRotationMasks.mat}  # initial transform
    DATA{Input/test_mask.nii.gz}  # fixedMask
    DATA{Input/rotation.test_mask.nii.gz}  # movingMask
    0.001  # (minimum?) StepLength
    Affine
    Mattes
    0.10  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    1x1x1  # shrinkFactors
    0.1x0.1x0.1 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_AffineRotationNoMasks
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
    none  # fixedMask
    none  # movingMask
    0.001  # (minimum?) StepLength
    Affine
    Mattes
    0.10  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    1x1x1  # shrinkFactors
    0.1x0.1x0.1 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_AffineTranslationMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_AffineTranslationMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineTranslationMasks.result.nii.gz
    --compareIntensityTolerance 7
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 777
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/translation.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineTranslationMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineTranslationMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_antsRegistrationTest_AffineTranslationMasks.mat}  # initial transform
    DATA{Input/test_mask.nii.gz}  # fixedMask
    DATA{Input/translation.test_mask.nii.gz}  # movingMask
    0.001  # (minimum?) StepLength
    Affine
    MeanSquares
    0.10  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    1x1x1  # shrinkFactors
    0.1x0.1x0.1 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_AffineTranslationNoMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_AffineTranslationNoMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineTranslationNoMasks.result.nii.gz
    --compareIntensityTolerance 7
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 777
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/translation.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineTranslationNoMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineTranslationNoMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_antsRegistrationTest_AffineTranslationMasks.mat}  # initial transform
    none  # fixedMask
    none  # movingMask
    0.001  # (minimum?) StepLength
    Affine
    MeanSquares
    0.10  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    1x1x1  # shrinkFactors
    0.1x0.1x0.1 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_RigidRotationHeadMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_RigidRotationHeadMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_RigidRotationHeadMasks.result.nii.gz
    --compareIntensityTolerance 7
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 777
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/rotation.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_RigidRotationHeadMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_RigidRotationHeadMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_RigidRotationHeadMasks.mat}  # initial transform
    none  # fixedMask
    none  # movingMask
    0.001  # (minimum?) StepLength
    Rigid
    Mattes
    0.20  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    1x1x1  # shrinkFactors
    0.1x0.1x0.1 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_RigidRotationNoMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_RigidRotationNoMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_RigidRotationNoMasks.result.nii.gz
    --compareIntensityTolerance 7
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 777
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/rotation.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_RigidRotationNoMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_RigidRotationNoMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_RigidRotationNoMasks.mat}  # initial transform
    none  # fixedMask
    none  # movingMask
    0.001  # (minimum?) StepLength
    Rigid
    Mattes
    0.20  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    1x1x1  # shrinkFactors
    0.1x0.1x0.1 # smoothingSigmas
  )
