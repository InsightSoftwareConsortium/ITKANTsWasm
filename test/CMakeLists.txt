itk_module_test()

set(ANTsWasmTests
  itkANTSRegistrationTest.cxx
  )

CreateTestDriver(ANTsWasm "${ANTsWasm-Test_LIBRARIES}" "${ANTsWasmTests}")

itk_add_test(NAME itkANTSRegistrationTest
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/itkANTSRegistrationTestOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkANTSRegistrationTestOutput.mha
  itkANTSRegistrationTest
    ${ITK_TEST_OUTPUT_DIR}/itkANTSRegistrationTestOutput.mha
  )
