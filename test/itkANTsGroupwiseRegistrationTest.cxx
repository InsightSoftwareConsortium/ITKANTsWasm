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

#include "itkANTsGroupwiseRegistration.h"

#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkTxtTransformIOFactory.h"
#include "itkTestingMacros.h"


int
itkANTsGroupwiseRegistrationTest(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cerr << "Missing parameters." << std::endl;
    std::cerr << "Usage: " << itkNameOfTestExecutableMacro(argv);
    std::cerr << " inputDirectory outputDirectory [numberOfFaces]";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  itk::TxtTransformIOFactory::RegisterOneFactory();

  std::string inDir = argv[1];
  std::string outDir = argv[2];

  unsigned numberOfFaces = 9;
  if (argc > 3)
  {
    numberOfFaces = std::stoul(argv[3]);
  }

  using ImageType = itk::Image<unsigned char, 2>;
  using FloatImageType = itk::Image<float, 2>;
  using FilterType = itk::ANTsGroupwiseRegistration<ImageType, FloatImageType, float>;
  typename FilterType::Pointer filter = FilterType::New();
  ITK_EXERCISE_BASIC_OBJECT_METHODS(filter, ANTsGroupwiseRegistration, ImageToImageFilter);

  itk::SimpleFilterWatcher        watcher(filter, "ANTs groupwise registration");
  std::vector<ImageType::Pointer> images;
  std::vector<ImageType *>        imagePointers;
  for (unsigned i = 1; i <= numberOfFaces; ++i)
  {
    std::string fileName = inDir + "/face" + std::to_string(i) + ".png";
    auto        image = itk::ReadImage<ImageType>(fileName);
    images.push_back(image);
    imagePointers.push_back(image.GetPointer());
  }

  filter->SetImageList(imagePointers);
  filter->DebugOn();
  filter->Update();

  auto templateImage = filter->GetTemplateImage();
  itk::WriteImage(templateImage, outDir + "/faceTemplate.nrrd");

  using TransformWriter = itk::TransformFileWriterTemplate<float>;
  TransformWriter::Pointer transformWriter = TransformWriter::New();
  for (unsigned i = 1; i < numberOfFaces; ++i)
  {
    auto fTransform = filter->GetTransform(i);
    transformWriter->SetFileName(outDir + "/face.tfm");
    transformWriter->SetInput(fTransform);
    transformWriter->Update();
  }

  return EXIT_SUCCESS;
}
