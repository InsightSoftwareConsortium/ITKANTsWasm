itk_module_test()

set(ANTsWasmTests
  itkANTSRegistrationTest.cxx
  itkANTSRegistrationBasicTests.cxx
  itkANTSGroupwiseRegistrationTestFaces.cxx
  itkANTSGroupwiseRegistrationTest3D.cxx
  )

CreateTestDriver(ANTsWasm "${ANTsWasm-Test_LIBRARIES}" "${ANTsWasmTests}")


itk_add_test(NAME itkANTSRegistrationBasicTests
  COMMAND ANTsWasmTestDriver
  itkANTSRegistrationBasicTests ${ITK_TEST_OUTPUT_DIR}
  )

itk_add_test(NAME itkANTsGroupwiseRegistrationTestFacesDirectSyN
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/itkANTsGroupwiseRegistrationTestFacesDirectSyN.png}
    ${ITK_TEST_OUTPUT_DIR}/itkANTsGroupwiseRegistrationTestFacesDirectSyN.nrrd
    --compareIntensityTolerance 10
    --compareRadiusTolerance 3
    --compareNumberOfPixelsTolerance 2000
  itkANTSGroupwiseRegistrationTestFaces
  DATA{${CMAKE_CURRENT_LIST_DIR}/Input/face10.png,REGEX:face[0-9]+\\.png}  # download data at build time
  ${ITK_TEST_OUTPUT_DIR}/itkANTsGroupwiseRegistrationTestFacesDirectSyN.nrrd
  10
  SyN
  )

# the rest are pairwise registration tests, mostly taken directly from reference ANTs (C++ and scripts)

itk_add_test(NAME antsRegistrationTest_AffineScaleMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_AffineScaleMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineScaleMasks.result.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/scale.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineScaleMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_AffineScaleMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_AffineScaleMasks.mat}  # initial transform
    DATA{Input/test_mask.nii.gz}  # fixedMask
    DATA{Input/scale.test_mask.nii.gz}  # movingMask
    0.001  # GradientStep
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
    0.001  # GradientStep
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
    0.001  # GradientStep
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
    0.001  # GradientStep
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
    0.001  # GradientStep
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
    0.001  # GradientStep
    Affine
    MeanSquares
    0.10  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    1x1x1  # shrinkFactors
    0.1x0.1x0.1 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_RigidAnisotropicMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_RigidAnisotropicMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_RigidAnisotropicMasks.result.nii.gz
    --compareIntensityTolerance 7
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 1500
  itkANTSRegistrationTest
    DATA{Input/ANON0006_20_T1_dbg_splayed.nii.gz}  # fixed image
    DATA{Input/ANON0006_20_T1_sag_twisted.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_RigidAnisotropicMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_RigidAnisotropicMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_antsRegistrationTest_RigidAnisotropicMasks.mat}  # initial transform
    none  # fixedMask
    none  # movingMask
    0.001  # GradientStep
    Rigid
    Mattes
    0.20  # samplingRate
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
    0.001  # GradientStep
    Rigid
    Mattes
    0.20  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    1x1x1  # shrinkFactors
    0.1x0.1x0.1 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_RigidRotationMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_RigidRotationMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_RigidRotationMasks.result.nii.gz
    --compareIntensityTolerance 7
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 777
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/rotation.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_RigidRotationMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_RigidRotationMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_RigidRotationHeadMasks.mat}  # initial transform
    DATA{Input/test_mask.nii.gz}  # fixedMask
    DATA{Input/rotation.test_mask.nii.gz}  # movingMask
    0.001  # GradientStep
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
    0.001  # GradientStep
    Rigid
    Mattes
    0.20  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    1x1x1  # shrinkFactors
    0.1x0.1x0.1 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_MSEAffineRotationMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_MSEAffineRotationMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_MSEAffineRotationMasks.result.nii.gz
    --compareIntensityTolerance 7
    --compareRadiusTolerance 2
    --compareNumberOfPixelsTolerance 1250
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/rotation.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_MSEAffineRotationMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_MSEAffineRotationMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_AffineRotationMasks.mat}  # initial transform
    DATA{Input/test_mask.nii.gz}  # fixedMask
    DATA{Input/rotation.test_mask.nii.gz}  # movingMask
    0.001  # GradientStep
    Affine
    MeanSquares
    0.10  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    1x1x1  # shrinkFactors
    0.1x0.1x0.1 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_SimilarityRotationMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_SimilarityRotationMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_SimilarityRotationMasks.result.nii.gz
    --compareIntensityTolerance 11
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/rotation.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_SimilarityRotationMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_SimilarityRotationMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_RigidRotationNoMasks.mat}  # initial transform
    DATA{Input/test_mask.nii.gz}  # fixedMask
    DATA{Input/rotation.test_mask.nii.gz}  # movingMask
    0.001  # GradientStep
    Similarity
    Mattes
    0.20  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    3x2x1  # shrinkFactors
    0x0x0 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_CCSimilarityRotationMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_CCSimilarityRotationMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_CCSimilarityRotationMasks.result.nii.gz
    --compareIntensityTolerance 11
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/rotation.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_CCSimilarityRotationMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_CCSimilarityRotationMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_RigidRotationNoMasks.mat}  # initial transform
    DATA{Input/test_mask.nii.gz}  # fixedMask
    DATA{Input/rotation.test_mask.nii.gz}  # movingMask
    0.001  # GradientStep
    Similarity
    CC
    0.20  # samplingRate
    5  # radius
    25x20x5  # affineIterations
    3x2x1  # shrinkFactors
    0x0x0 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_MSESimilarityRotationMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_MSESimilarityRotationMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_MSESimilarityRotationMasks.result.nii.gz
    --compareIntensityTolerance 11
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/rotation.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_MSESimilarityRotationMasks.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_MSESimilarityRotationMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_RigidRotationNoMasks.mat}  # initial transform
    DATA{Input/test_mask.nii.gz}  # fixedMask
    DATA{Input/rotation.test_mask.nii.gz}  # movingMask
    0.001  # GradientStep
    Similarity
    MeanSquares
    0.20  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    3x2x1  # shrinkFactors
    0x0x0 # smoothingSigmas
  )

itk_add_test(NAME antsRegistrationTest_MSESimilarityRotationMasks_Float
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_MSESimilarityRotationMasks_Float.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_MSESimilarityRotationMasks_Float.result.nii.gz
    --compareIntensityTolerance 11
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/rotation.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_MSESimilarityRotationMasks_Float.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_MSESimilarityRotationMasks_Float.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_RigidRotationNoMasks_Float.mat}  # initial transform
    DATA{Input/test_mask.nii.gz}  # fixedMask
    DATA{Input/rotation.test_mask.nii.gz}  # movingMask
    0.001  # GradientStep
    Similarity
    MeanSquares
    0.20  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    3x2x1  # shrinkFactors
    0x0x0 # smoothingSigmas
    0  # randomSeed (0 means do not set)
    Mattes  # synMetric
    100x70x20  # synIterations  
    --float    
  )

itk_add_test(NAME antsRegistrationTest_SyNScaleNoMasks
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_SyNScaleNoMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_SyNScaleNoMasks.result.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/scale.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_SyNScaleNoMasks.h5  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_SyNScaleNoMasks.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_AffineScaleMasks.mat}  # initial transform
    none  # fixedMask
    none  # movingMask
    0.25  # GradientStep
    SyNOnly
    Mattes  # affineMetric
    0.20  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    3x2x1  # shrinkFactors
    2x1x0 # smoothingSigmas
    0  # randomSeed (0 means do not set)
    Mattes  # synMetric
    100x70x20  # synIterations
  )

itk_add_test(NAME antsRegistrationTest_SyNScaleNoMasks_Float
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_SyNScaleNoMasks_Float.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_SyNScaleNoMasks_Float.result.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/scale.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_SyNScaleNoMasks_Float.h5  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_SyNScaleNoMasks_Float.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_AffineScaleMasks_Float.mat}  # initial transform
    none  # fixedMask
    none  # movingMask
    0.25  # GradientStep
    SyNOnly
    Mattes  # affineMetric
    0.20  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    3x2x1  # shrinkFactors
    2x1x0 # smoothingSigmas
    0  # randomSeed (0 means do not set)
    Mattes  # synMetric
    100x70x20  # synIterations
    --float
  )

itk_add_test(NAME antsRegistrationTest_SyNScaleNoMasks_2stage
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_SyNScaleNoMasks_Float.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_SyNScaleNoMasks_2stage.result.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/scale.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_SyNScaleNoMasks_2stage.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistrationTest_SyNScaleNoMasks_2stage.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_AffineScaleMasks_Float.mat}  # initial transform
    none  # fixedMask
    none  # movingMask
    0.25  # GradientStep
    SyN
    Mattes  # affineMetric
    0.20  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    3x2x1  # shrinkFactors
    2x1x0 # smoothingSigmas
    0  # randomSeed (0 means do not set)
    Mattes  # synMetric
    100x70x20  # synIterations
    --float
    0  # collapseTransforms
    1  # maskAllStages
  )

itk_add_test(NAME antsRegistration_SyNCC
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_SyNScaleNoMasks_Float.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistration_SyNCC.result.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/scale.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistration_SyNCC.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistration_SyNCC.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_AffineScaleMasks_Float.mat}  # initial transform
    none  # fixedMask
    none  # movingMask
    0.25  # GradientStep
    SyNCC
    Irrelevant  # affineMetric
    0.20  # samplingRate
    4  # ccRadius
    25x20x5  # affineIterations
    3x2x1  # shrinkFactors
    2x1x0 # smoothingSigmas
    0  # randomSeed (0 means do not set)
    Irrelevant  # synMetric
    100x70x20  # synIterations
    --float
    0  # collapseTransforms
    1  # maskAllStages
  )

itk_add_test(NAME antsRegistration_SyNRA
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/antsRegistrationTest_SyNScaleNoMasks_Float.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/antsRegistration_SyNRA.result.nii.gz
    --compareIntensityTolerance 9
    --compareRadiusTolerance 1
    --compareNumberOfPixelsTolerance 1000
  itkANTSRegistrationTest
    DATA{Input/test.nii.gz}  # fixed image
    DATA{Input/scale.test.nii.gz}  # moving image
    ${ITK_TEST_OUTPUT_DIR}/antsRegistration_SyNRA.tfm  # output transform
    ${ITK_TEST_OUTPUT_DIR}/antsRegistration_SyNRA.result.nii.gz  # moving image warped to fixed space
    DATA{Input/Initializer_0.05_antsRegistrationTest_AffineScaleMasks_Float.mat}  # initial transform
    none  # fixedMask
    none  # movingMask
    0.25  # GradientStep
    SyNRA
    Mattes  # affineMetric
    0.20  # samplingRate
    200  # numberOfBins
    25x20x5  # affineIterations
    3x2x1  # shrinkFactors
    2x1x0 # smoothingSigmas
    0  # randomSeed (0 means do not set)
    Mattes  # synMetric
    100x70x20  # synIterations
    --float
    0  # collapseTransforms
    1  # maskAllStages
  )
