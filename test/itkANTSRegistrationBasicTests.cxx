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
#include "itkSimpleFilterWatcher.h"
#include "itkImageRegionIterator.h"
#include "itkTestingMacros.h"

namespace
{
template <typename TImage>
void setRegionToValue(TImage * image, const typename TImage::RegionType region, typename TImage::PixelType value)
{
  itk::ImageRegionIterator<TImage> imageIterator(image, region);
  while (!imageIterator.IsAtEnd())
  {
    imageIterator.Set(value);
    ++imageIterator;
  }
};

template <typename FixedImageType, typename MovingImageType>
int
testFilter()
{
  using FilterType = itk::ANTSRegistration<FixedImageType, MovingImageType>;
  typename FilterType::Pointer filter = FilterType::New();
  ITK_EXERCISE_BASIC_OBJECT_METHODS(filter, ANTSRegistration, ProcessObject);

  // Create input images to avoid test dependencies.
  typename FixedImageType::SizeType size;
  size.Fill(128);
  typename FixedImageType::Pointer fixedImage = FixedImageType::New();
  fixedImage->SetRegions(size);
  fixedImage->Allocate();
  fixedImage->FillBuffer(11.1f);
  typename FixedImageType::RegionType region = fixedImage->GetLargestPossibleRegion();
  region.ShrinkByRadius(20);
  setRegionToValue(fixedImage.GetPointer(), region, 22.2f);

  typename MovingImageType::Pointer movingImage = MovingImageType::New();
  movingImage->SetRegions(size);
  movingImage->Allocate();
  movingImage->FillBuffer(44);
  region.SetIndex(0, region.GetIndex(0) + 10); // shift the rectangle
  setRegionToValue(movingImage.GetPointer(), region, 55);

  itk::SimpleFilterWatcher watcher(filter, "ANTs registration");

  filter->SetFixedImage(fixedImage);
  filter->SetMovingImage(movingImage);
  filter->SetTypeOfTransform("Affine");
  filter->Update();

  auto filterOutput = filter->GetForwardTransform();
  std::cout << *filterOutput << std::endl;
  typename MovingImageType::Pointer movingResampled = filter->GetWarpedMovingImage();
  std::cout << *movingResampled << std::endl;
  auto inverseTransform = filter->GetInverseTransform(); // Affine should be invertible
  std::cout << *inverseTransform << std::endl;
  typename FixedImageType::Pointer fixedResampled = filter->GetWarpedFixedImage();
  std::cout << *fixedResampled << std::endl;

  return EXIT_SUCCESS;
}
} // namespace


int
itkANTSRegistrationBasicTests(int, char ** const)
{
  std::cout << "\nTesting: fixed and moving image are of the same pixel type" << std::endl;
  using Float2DImageType = itk::Image<float, 2>;
  int retVal = testFilter<Float2DImageType, Float2DImageType>();
  if (retVal != EXIT_SUCCESS)
  {
    return retVal;
  }

  std::cout << "\nTesting: fixed and moving image are of different pixel type" << std::endl;
  using Short2DImageType = itk::Image<short, 2>;
  return testFilter<Float2DImageType, Short2DImageType>();
}
