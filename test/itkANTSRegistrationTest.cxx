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

#include "itkCommand.h"
#include "itkImageFileWriter.h"
#include "itkImageDuplicator.h"
#include "itkTestingMacros.h"

namespace
{
class ShowProgress : public itk::Command
{
public:
  itkNewMacro(ShowProgress);

  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    Execute((const itk::Object *)caller, event);
  }

  void
  Execute(const itk::Object * caller, const itk::EventObject & event) override
  {
    if (!itk::ProgressEvent().CheckEvent(&event))
    {
      return;
    }
    const auto * processObject = dynamic_cast<const itk::ProcessObject *>(caller);
    if (!processObject)
    {
      return;
    }
    std::cout << " " << processObject->GetProgress();
  }
};
} // namespace


constexpr unsigned int Dimension = 2;
using PixelType = float;
using ImageType = itk::Image<PixelType, Dimension>;

auto setRegionToValue = [](ImageType::Pointer image, ImageType::RegionType region, PixelType value) {
  itk::ImageRegionIterator<ImageType> imageIterator(image, region);
  while (!imageIterator.IsAtEnd())
  {
    imageIterator.Set(value);
    ++imageIterator;
  }
};

int
itkANTSRegistrationTest(int argc, char * argv[])
{
  if (argc < 2)
  {
    std::cerr << "Missing parameters." << std::endl;
    std::cerr << "Usage: " << itkNameOfTestExecutableMacro(argv);
    std::cerr << " outputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }
  const char * outputImageFileName = argv[1];


  using FilterType = itk::ANTSRegistration<ImageType, ImageType>;
  FilterType::Pointer filter = FilterType::New();

  ITK_EXERCISE_BASIC_OBJECT_METHODS(filter, ANTSRegistration, ProcessObject);

  // Create input image to avoid test dependencies.
  ImageType::SizeType size;
  size.Fill(128);
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(size);
  image->Allocate();
  image->FillBuffer(1.1f);
  ImageType::RegionType region = image->GetLargestPossibleRegion();
  region.ShrinkByRadius(20);
  setRegionToValue(image, region, 0.2f);

  using DuplicatorType = itk::ImageDuplicator<ImageType>;
  auto duplicator = DuplicatorType::New();
  duplicator->SetInputImage(image);
  duplicator->Update();
  ImageType::Pointer clonedImage = duplicator->GetOutput();
  region.SetIndex(0, region.GetIndex(0) + 10); // shift the rectangle
  clonedImage->FillBuffer(0.3f);
  setRegionToValue(clonedImage, region, 1.2f);

  ShowProgress::Pointer showProgress = ShowProgress::New();
  filter->AddObserver(itk::ProgressEvent(), showProgress);
  filter->SetFixedImage(image);
  filter->SetMovingImage(clonedImage);
  filter->SetTypeOfTransform("Affine");
  filter->Update();
  auto filterOutput = filter->GetForwardTransform();
  std::cout << "\nForwardTransform: " << filterOutput << std::endl;

  using WriterType = itk::ImageFileWriter<ImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputImageFileName);
  writer->SetInput(image);
  writer->SetUseCompression(true);

  ITK_TRY_EXPECT_NO_EXCEPTION(writer->Update());

  std::cout << "Test finished." << std::endl;
  return EXIT_SUCCESS;
}
