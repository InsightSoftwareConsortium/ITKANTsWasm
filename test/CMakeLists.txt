itk_module_test()

set(ANTsWasmTests
  itkMinimalStandardRandomVariateGeneratorTest.cxx
  itkANTSRegistrationTest.cxx
  )

CreateTestDriver(ANTsWasm "${ANTsWasm-Test_LIBRARIES}" "${ANTsWasmTests}")

itk_add_test(NAME itkMinimalStandardRandomVariateGeneratorTest
  COMMAND ANTsWasmTestDriver itkMinimalStandardRandomVariateGeneratorTest
  )

itk_add_test(NAME itkANTSRegistrationTest
  COMMAND ANTsWasmTestDriver
    --compare
    DATA{Baseline/itkANTSRegistrationTestOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkANTSRegistrationTestOutput.mha
  itkANTSRegistrationTest
    ${ITK_TEST_OUTPUT_DIR}/itkANTSRegistrationTestOutput.mha
  )
