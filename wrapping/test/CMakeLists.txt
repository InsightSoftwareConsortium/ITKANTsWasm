set(test_input_dir ${${itk-module}_SOURCE_DIR}/test/Input)
set(test_baseline_dir ${${itk-module}_SOURCE_DIR}/test/Baseline)

itk_python_expression_add_test(NAME PythonInstantiateANTSRegistrationFloat2D
  EXPRESSION "reg = itk.ANTSRegistration[itk.Image[itk.F, 2], itk.Image[itk.F, 2]].New()"
  )

itk_python_add_test(NAME PythonANTSRegistrationTest
  TEST_DRIVER_ARGS
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777
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
