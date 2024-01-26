/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkANTSRegistration.h"

#include "itkImageFileWriter.h"
#include "itkMatlabTransformIOFactory.h"
#include "itkTxtTransformIOFactory.h"
#include "itkTestingMacros.h"

namespace
{
template <unsigned Dimension>
int
doTest(int argc, char * argv[])
{
  const char * fixedImageFileName = argv[1];
  const char * movingImageFileName = argv[2];
  const char * outTransformFileName = argv[3];

  itk::MatlabTransformIOFactory::RegisterOneFactory();
  itk::TxtTransformIOFactory::RegisterOneFactory();

  using ImageType = itk::Image<float, Dimension>;
  using LabelImageType = itk::Image<unsigned char, Dimension>;
  using FilterType = itk::ANTSRegistration<ImageType, ImageType>;
  typename FilterType::Pointer filter = FilterType::New();

  typename ImageType::Pointer fixedImage;
  ITK_TRY_EXPECT_NO_EXCEPTION(fixedImage = itk::ReadImage<ImageType>(fixedImageFileName));

  typename ImageType::Pointer movingImage;
  ITK_TRY_EXPECT_NO_EXCEPTION(movingImage = itk::ReadImage<ImageType>(movingImageFileName));

  if (argc > 5)
  {
    itk::TransformFileReader::Pointer transformReader = itk::TransformFileReader::New();
    transformReader->SetFileName(argv[5]);
    ITK_TRY_EXPECT_NO_EXCEPTION(transformReader->Update());
    const itk::TransformFileReader::TransformListType * initialTransformList = transformReader->GetTransformList();
    ITK_TEST_EXPECT_EQUAL(initialTransformList->size(), 1);

    auto firstTransform = initialTransformList->begin();
    bool transformSupported = false;
    if constexpr (Dimension == 3)
    {
      if (!strcmp((*firstTransform)->GetNameOfClass(), "VersorRigid3DTransform"))
      {
        using VersorType = itk::VersorRigid3DTransform<double>;
        typename VersorType::Pointer versor = static_cast<VersorType *>((*firstTransform).GetPointer());
        versor->Print(std::cout);
        filter->SetInitialTransform(versor);
        transformSupported = true;
      }
    }
    if (!strcmp((*firstTransform)->GetNameOfClass(), "CompositeTransform"))
    {
      using CompositeType = itk::CompositeTransform<double, Dimension>;
      typename CompositeType::Pointer composite = static_cast<CompositeType *>((*firstTransform).GetPointer());
      composite->Print(std::cout);
      filter->SetInitialTransform(composite);
    }
    else if (!strcmp((*firstTransform)->GetNameOfClass(), "AffineTransform"))
    {
      using AffineType = itk::AffineTransform<double, Dimension>;
      typename AffineType::Pointer affine = static_cast<AffineType *>((*firstTransform).GetPointer());
      affine->Print(std::cout);
      filter->SetInitialTransform(affine);
    }
    else if (!transformSupported)
    {
      std::cout << "Unsupported initial transform type: " << (*firstTransform)->GetNameOfClass() << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (argc > 6 && std::string(argv[6]) != "none")
  {
    typename LabelImageType::Pointer fixedMask;
    ITK_TRY_EXPECT_NO_EXCEPTION(fixedMask = itk::ReadImage<LabelImageType>(argv[6]));
    filter->SetFixedMask(fixedMask);
  }
  if (argc > 7 && std::string(argv[7]) != "none")
  {
    typename LabelImageType::Pointer movingMask;
    ITK_TRY_EXPECT_NO_EXCEPTION(movingMask = itk::ReadImage<LabelImageType>(argv[7]));
    filter->SetMovingMask(movingMask);
  }
  if (argc > 8)
  {
    double gradientStep = std::stod(argv[8]);
    filter->SetGradientStep(gradientStep);
  }

  filter->SetFixedImage(fixedImage);
  filter->SetMovingImage(movingImage);
  filter->SetTypeOfTransform("Affine");
  ITK_TRY_EXPECT_NO_EXCEPTION(filter->Update());

  // debug
  auto filterOutput = filter->GetForwardTransform();
  std::cout << "\nForwardTransform: " << *filterOutput << std::endl;

  itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
  transformWriter->SetInput(filterOutput);
  transformWriter->SetFileName(outTransformFileName);
  ITK_TRY_EXPECT_NO_EXCEPTION(transformWriter->Update());

  if (argc > 4)
  {
    ITK_TRY_EXPECT_NO_EXCEPTION(itk::WriteImage(filter->GetWarpedMovingImage(), argv[4]));
  }

  std::cout << "Test finished." << std::endl;
  return EXIT_SUCCESS;
}

} // namespace


int
itkANTSRegistrationTest(int argc, char * argv[])
{
  if (argc < 4)
  {
    std::cerr << "Missing parameters." << std::endl;
    std::cerr << "Usage: " << itkNameOfTestExecutableMacro(argv);
    std::cerr << " fixedImage movingImage outTransform";
    std::cerr << " [outMovingResampledToFixedSpace] [initialTransform]";
    std::cerr << " [fixedImageMask] [movingImageMask]";
    std::cerr << " [gradientStep]";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  // determine dimension
  using ImageType = itk::Image<unsigned char, 2>;
  auto imageReader = itk::ImageFileReader<ImageType>::New();
  imageReader->SetFileName(argv[1]);
  ITK_TRY_EXPECT_NO_EXCEPTION(imageReader->UpdateOutputInformation());
  unsigned dimension = imageReader->GetImageIO()->GetNumberOfDimensions();

  switch (dimension)
  {
    case 2:
      return doTest<2>(argc, argv);
    case 3:
      return doTest<3>(argc, argv);
    case 4:
      return doTest<4>(argc, argv);
    default:
      std::cerr << "Unsupported image dimension: " << dimension;
      return EXIT_FAILURE;
  }
}
