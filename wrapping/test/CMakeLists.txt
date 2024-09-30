set(test_input_dir ${${itk-module}_SOURCE_DIR}/test/Input)
set(test_baseline_dir ${${itk-module}_SOURCE_DIR}/test/Baseline)

itk_python_expression_add_test(NAME PythonInstantiateANTSRegistrationFloat2D
  EXPRESSION "reg = itk.ANTSRegistration[itk.Image[itk.F, 2], itk.Image[itk.F, 2], itk.D].New()"
  )

itk_python_add_test(NAME PythonANTSRegistrationTest
  TEST_DRIVER_ARGS
  --compareIntensityTolerance 11
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1000
  --compare
    DATA{${test_baseline_dir}/antsRegistrationTest_AffineTranslationNoMasks.result.nii.gz}
    ${ITK_TEST_OUTPUT_DIR}/PythonANTSRegistrationTest_AffineTranslationNoMasks.result.nii.gz
  COMMAND PythonANTSRegistrationTest.py
    --fixed-image DATA{${test_input_dir}/test.nii.gz}
    --moving-image DATA{${test_input_dir}/translation.test.nii.gz}
    --initial-transform DATA{${test_input_dir}/Initializer_antsRegistrationTest_AffineTranslationMasks.mat}
    --output-transform ${ITK_TEST_OUTPUT_DIR}/PythonANTSRegistrationTest_AffineTranslationNoMasks.tfm
    --resampled-moving ${ITK_TEST_OUTPUT_DIR}/PythonANTSRegistrationTest_AffineTranslationNoMasks.result.nii.gz
  )

itk_python_expression_add_test(NAME PythonInstantiateANTSGroupwiseRegistrationShort3Dfloat
  EXPRESSION "reg = itk.ANTSGroupwiseRegistration[itk.Image[itk.SS, 3], itk.Image[itk.F, 3], itk.F].New()"
  )

itk_python_add_test(NAME PythonANTsGroupwiseFaces
  TEST_DRIVER_ARGS
  --compare
    DATA{${test_baseline_dir}/PythonANTsGroupwiseFaces.png}
    ${ITK_TEST_OUTPUT_DIR}/PythonANTsGroupwiseFaces.nrrd
    --compareIntensityTolerance 10
    --compareRadiusTolerance 3
    --compareNumberOfPixelsTolerance 2000
  COMMAND PythonANTsGroupwiseFaces.py
    -o ${ITK_TEST_OUTPUT_DIR}/PythonANTsGroupwiseFaces.nrrd
    -i DATA{${test_input_dir}/face1.png}
    -i DATA{${test_input_dir}/face2.png}
    -i DATA{${test_input_dir}/face3.png}
    -i DATA{${test_input_dir}/face4.png}
    -i DATA{${test_input_dir}/face5.png}
    -i DATA{${test_input_dir}/face6.png}
    -i DATA{${test_input_dir}/face7.png}
    -i DATA{${test_input_dir}/face8.png}
    -i DATA{${test_input_dir}/face9.png}
    -i DATA{${test_input_dir}/face10.png}
  )

itk_python_add_test(NAME PythonANTsGroupwiseBrainSlices
  TEST_DRIVER_ARGS
  --compare
    DATA{${test_baseline_dir}/PythonANTsGroupwiseBrainSlices.nrrd}
    ${ITK_TEST_OUTPUT_DIR}/PythonANTsGroupwiseBrainSlices.nrrd
    --compareIntensityTolerance 10
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 2000
  COMMAND PythonANTsGroupwiseBrainSlices.py
    -o ${ITK_TEST_OUTPUT_DIR}/PythonANTsGroupwiseBrainSlices.nrrd
    -i DATA{${test_input_dir}/OASIS-TRT-20-10Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-11Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-12Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-13Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-14Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-15Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-16Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-17Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-18Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-19Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-20Slice121.nii.gz}
  )

itk_python_add_test(NAME PythonANTsGroupwisePaths
  TEST_DRIVER_ARGS
  --compare
    DATA{${test_baseline_dir}/PythonANTsGroupwiseBrainSlices.nrrd}
    ${ITK_TEST_OUTPUT_DIR}/PythonANTsGroupwisePaths.nrrd
    --compareIntensityTolerance 10
    --compareRadiusTolerance 5
    --compareNumberOfPixelsTolerance 2000
  COMMAND PythonANTsGroupwisePaths.py
    -o ${ITK_TEST_OUTPUT_DIR}/PythonANTsGroupwisePaths.nrrd
    -i DATA{${test_input_dir}/OASIS-TRT-20-10Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-11Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-12Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-13Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-14Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-15Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-16Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-17Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-18Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-19Slice121.nii.gz}
    -i DATA{${test_input_dir}/OASIS-TRT-20-20Slice121.nii.gz}
  )
