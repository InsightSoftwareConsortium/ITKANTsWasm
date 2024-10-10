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
#ifndef itkANTSGroupwiseBuildTemplate_hxx
#define itkANTSGroupwiseBuildTemplate_hxx

#include <sstream>

#include "itkPrintHelper.h"
#include "itkANTSGroupwiseBuildTemplate.h" // needed by Visual Studio for code completion
#include "itkImageDuplicator.h"
#include "itkWeightedAddImageFilter.h"
#include "itkAverageAffineTransformFunction.h"
#include "itkLaplacianSharpeningImageFilter.h"

namespace // anonymous namespace for debug functions
{
template <typename TransformType>
void
WriteTransform(const TransformType * transform, std::string filename)
{
  constexpr unsigned Dimension = TransformType::OutputSpaceDimension;
  using Affine3D = itk::AffineTransform<typename TransformType::ParametersValueType, 3>;
  using Affine2D = itk::AffineTransform<typename TransformType::ParametersValueType, 2>;
  using TransformWriterType = itk::TransformFileWriterTemplate<typename TransformType::ParametersValueType>;
  typename TransformWriterType::Pointer tWriter = TransformWriterType::New();
  tWriter->SetFileName(filename);

  if (Dimension == 2 && std::string(transform->GetNameOfClass()) == "AffineTransform")
  {
    typename Affine2D::ConstPointer t2d = dynamic_cast<const Affine2D *>(transform);
    // convert into affine which Slicer can read
    typename Affine3D::MatrixType      m;
    typename Affine3D::TranslationType t;
    t.Fill(0);
    typename Affine3D::InputPointType c;
    c.Fill(0);
    for (unsigned int i = 0; i < Dimension; i++)
    {
      for (unsigned int j = 0; j < Dimension; j++)
      {
        m[i][j] = t2d->GetMatrix()(i, j);
      }
      t[i] = t2d->GetMatrix()(i, Dimension);
      c[i] = t2d->GetCenter()[i];
    }
    m(Dimension, Dimension) = 1.0;

    typename Affine3D::Pointer aTr = Affine3D::New();
    aTr->SetCenter(c);
    aTr->SetMatrix(m);
    aTr->SetOffset(t);
    tWriter->SetInput(aTr);
  }
  else
  {
    tWriter->SetInput(transform);
  }
  tWriter->Update();
}
} // namespace

namespace itk
{

template <typename TImage, typename TTemplateImage, typename TParametersValueType>
ANTSGroupwiseBuildTemplate<TImage, TTemplateImage, TParametersValueType>::ANTSGroupwiseBuildTemplate()
{
  this->SetPrimaryInputName("InitialTemplate");
  this->SetPrimaryOutputName("OptimizedImage");

  this->SetInput(0, TemplateImageType::New()); // empty initial template

  this->GetMultiThreader()->SetMaximumNumberOfThreads(1); // registrations are already multi-threaded
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
void
ANTSGroupwiseBuildTemplate<TImage, TTemplateImage, TParametersValueType>::PrintSelf(std::ostream & os,
                                                                                   Indent         indent) const
{
  using namespace print_helper;
  Superclass::PrintSelf(os, indent);

  os << indent << "GradientStep: " << this->m_GradientStep << '\n';
  os << indent << "BlendingWeight: " << this->m_BlendingWeight << '\n';
  os << indent << "UseNoRigid: " << (this->m_UseNoRigid ? "On" : "Off") << '\n';
  os << indent << "Iterations: " << this->m_Iterations << '\n';
  os << indent << "Weights: " << this->m_Weights << '\n';
  os << indent << "PathList: " << this->m_PathList << '\n';

  os << indent << "ImageList: " << '\n';
  unsigned i = 0;
  for (const auto & image : this->m_ImageList)
  {
    os << indent.GetNextIndent() << "Image" << i++ << ": ";
    // image->Print(os, indent.GetNextIndent());
    os << image.GetPointer() << '\n';
  }

  os << indent << "PairwiseRegistration: ";
  if (this->m_PairwiseRegistration)
  {
    this->m_PairwiseRegistration->Print(os, indent.GetNextIndent());
  }
  else
  {
    os << "nullptr" << std::endl;
  }
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
void
ANTSGroupwiseBuildTemplate<TImage, TTemplateImage, TParametersValueType>::VerifyInputInformation() const
{
  if (m_ImageList.empty() && m_PathList.empty())
  {
    itkExceptionMacro("Neither input images nor filename paths are provided.");
  }

  if (m_ImageList.empty() == m_PathList.empty())
  {
    itkExceptionMacro("Either input images OR filename paths must be provided, but not both.");
  }

  std::size_t imageListSize = std::max(m_ImageList.size(), m_PathList.size());

  if (m_Weights.size() > 0 && m_Weights.size() != imageListSize)
  {
    itkExceptionMacro("The number of weights is different from the number of images.");
  }

  if (imageListSize < 2)
  {
    itkExceptionMacro("At least two input images are required.");
  }
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
DataObject::Pointer
ANTSGroupwiseBuildTemplate<TImage, TTemplateImage, TParametersValueType>::MakeOutput(DataObjectPointerArraySizeType)
{
  typename OutputTransformType::Pointer ptr;
  Self::MakeOutputTransform(ptr);
  typename DecoratedOutputTransformType::Pointer decoratedOutputTransform = DecoratedOutputTransformType::New();
  decoratedOutputTransform->Set(ptr);
  return decoratedOutputTransform;
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
auto
ANTSGroupwiseBuildTemplate<TImage, TTemplateImage, TParametersValueType>::DuplicateImage(const TemplateImageType * image)
  -> typename TemplateImageType::Pointer
{
  using DuplicatorType = ImageDuplicator<TemplateImageType>;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(image);
  duplicator->Update();
  return duplicator->GetOutput();
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
template <typename TOutputImage, typename TInputImage>
auto
itk::ANTSGroupwiseBuildTemplate<TImage, TTemplateImage, TParametersValueType>::ResampleToTarget(
  const TInputImage *                  input,
  const TemplateImageType *            target,
  typename TransformType::ConstPointer transform) -> typename TOutputImage::Pointer
{
  using ResampleFilterType = ResampleImageFilter<TInputImage, TOutputImage, ParametersValueType, ParametersValueType>;
  typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
  resampleFilter->SetInput(input);
  if (transform)
  {
    resampleFilter->SetTransform(transform);
  }
  resampleFilter->SetOutputParametersFromImage(target);
  resampleFilter->Update();
  typename TOutputImage::Pointer result = resampleFilter->GetOutput();
  result->DisconnectPipeline();
  return result;
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
template <typename TTempImage>
typename TTempImage::Pointer
ANTSGroupwiseBuildTemplate<TImage, TTemplateImage, TParametersValueType>::ScaleAndAdd(
  typename TTempImage::Pointer   temp,
  const TTempImage *             image,
  typename TTempImage::PixelType weight)
{
  using MulType = MultiplyImageFilter<TTempImage, TTempImage, TTempImage>;
  typename MulType::Pointer mul = MulType::New();
  mul->SetInput1(image);
  mul->SetConstant2(weight);
  mul->Update();

  using AddImageFilterType = AddImageFilter<TTempImage, TTempImage, TTempImage>;
  typename AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
  addImageFilter->SetInput1(temp);
  addImageFilter->SetInput2(mul->GetOutput());
  addImageFilter->SetInPlace(true);
  addImageFilter->Update();
  temp = addImageFilter->GetOutput();
  temp->DisconnectPipeline();

  return temp;
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
void
ANTSGroupwiseBuildTemplate<TImage, TTemplateImage, TParametersValueType>::GenerateOutputInformation()
{
  if (!m_PairwiseRegistration) // a custom pairwise registration is not set
  {
    m_PairwiseRegistration = PairwiseType::New();
    m_PairwiseRegistration->SetTypeOfTransform("SyN");
    // m_PairwiseRegistration->DebugOn();
  }

  std::size_t imageListSize = std::max(m_ImageList.size(), m_PathList.size());
  if (m_Weights.empty())
  {
    m_Weights.resize(imageListSize, 1.0 / imageListSize);
  }
  else // normalize to sum to 1
  {
    m_Weights.resize(imageListSize, 1.0);
    TParametersValueType sum = 0;
    for (const auto & weight : m_Weights)
    {
      sum += weight;
    }
    for (auto & weight : m_Weights)
    {
      weight /= sum;
    }
  }
  m_TransformList.resize(imageListSize, nullptr);

  typename TemplateImageType::Pointer initialTemplate = dynamic_cast<TemplateImageType *>(this->GetInput(0));
  if (!initialTemplate)
  {
    itkExceptionMacro("Initial template must be a float-pixel image.");
  }

  TTemplateImage * outputPtr = this->GetOutput();
  if (initialTemplate->GetLargestPossibleRegion().GetNumberOfPixels() > 0)
  {
    // copy the initial template into the output (the current average template)
    outputPtr->CopyInformation(initialTemplate);
    outputPtr->SetLargestPossibleRegion(initialTemplate->GetLargestPossibleRegion());
  }
  else
  {
    // get the output metadata from the first image
    if (m_PathList.empty())
    {
      outputPtr->CopyInformation(m_ImageList[0]);
      outputPtr->SetLargestPossibleRegion(m_ImageList[0]->GetLargestPossibleRegion());
    }
    else
    {
      if (m_KeepTransforms)
      {
        itkExceptionMacro("KeepTransforms defeats the memory-saving purpose of reading images from files.");
      }
      using ReaderType = ImageFileReader<TImage>;
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName(m_PathList[0]);
      reader->UpdateOutputInformation();
      outputPtr->CopyInformation(reader->GetOutput());
      outputPtr->SetLargestPossibleRegion(reader->GetOutput()->GetLargestPossibleRegion());
    }
  }
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
void
ANTSGroupwiseBuildTemplate<TImage, TTemplateImage, TParametersValueType>::GenerateData()
{
  this->UpdateProgress(0.0);

  std::size_t imageListSize = std::max(m_ImageList.size(), m_PathList.size());
  this->GenerateOutputInformation();
  TTemplateImage * outputPtr = this->GetOutput();

  typename TemplateImageType::Pointer initialTemplate = dynamic_cast<TemplateImageType *>(this->GetInput(0));
  if (initialTemplate->GetLargestPossibleRegion().GetNumberOfPixels() > 0)
  {
    // copy the initial template into the output (the current average template)
    auto region = initialTemplate->GetLargestPossibleRegion();
    outputPtr->SetRegions(region);
    outputPtr->Allocate(false);
    ImageAlgorithm::Copy(initialTemplate.GetPointer(), outputPtr, region, region);
  }
  else
  {
    // create a simple average of all the input images
    typename TemplateImageType::Pointer average = TemplateImageType::New();
    average->CopyInformation(outputPtr); // GenerateOutputInformation() has already set the output metadata
    average->SetRegions(outputPtr->GetLargestPossibleRegion());
    average->Allocate(true); // initialize to zero

    for (unsigned i = 0; i < imageListSize; ++i)
    {
      typename ImageType::Pointer image;
      if (m_PathList.empty())
      {
        image = m_ImageList[i];
      }
      else
      {
        image = ReadImage<TImage>(m_PathList[i]);
      }
      typename TemplateImageType::Pointer resampledImage =
        ResampleToTarget<TemplateImageType, ImageType>(image, average, nullptr);
      average = ScaleAndAdd<TemplateImageType>(average, resampledImage, m_Weights[i]);
    } // both image and resampledImage should go out of scope here

    outputPtr->SetRegions(average->GetLargestPossibleRegion());
    outputPtr->SetPixelContainer(average->GetPixelContainer());
    outputPtr->Modified();
  }
  this->UpdateProgress(0.01);
  // WriteImage(this->GetOutput(), "initialTemplate.nrrd"); // debug

  float progressStep = 0.98 / (m_Iterations * imageListSize);

  for (unsigned i = 0; i < m_Iterations; ++i)
  {
    using AffineType = AffineTransform<ParametersValueType, ImageDimension>;
    using DisplacementTransformType = DisplacementFieldTransform<ParametersValueType, ImageDimension>;
    using DisplacementImageType = typename DisplacementTransformType::DisplacementFieldType;

    typename TemplateImageType::Pointer xavg = this->DuplicateImage(this->GetOutput());
    // WriteImage(xavg, "xavg" + std::to_string(i) + ".nrrd"); // debug
    typename TemplateImageType::Pointer xavgNew = TemplateImageType::New();
    xavgNew->CopyInformation(xavg);
    xavgNew->SetRegions(xavg->GetLargestPossibleRegion());
    xavgNew->Allocate(true); // initialize to zero
    typename DisplacementImageType::Pointer wavg{ nullptr };

    std::vector<typename AffineType::ConstPointer> affineList(imageListSize, nullptr);
    for (unsigned k = 0; k < imageListSize; ++k)
    {
      typename ImageType::Pointer image;
      if (m_PathList.empty())
      {
        image = m_ImageList[k];
      }
      else
      {
        image = ReadImage<TImage>(m_PathList[k]);
      }

      m_PairwiseRegistration->SetFixedImage(xavg);
      m_PairwiseRegistration->SetMovingImage(image);
      m_PairwiseRegistration->Update();

      const CompositeTransformType * compositeTransform = m_PairwiseRegistration->GetForwardTransform();
      // WriteTransform(compositeTransform, "tc" + std::to_string(i) + "_" + std::to_string(k) + ".tfm"); // debug
      affineList[k] = dynamic_cast<const AffineType *>(compositeTransform->GetFrontTransform());
      auto dfTransform = dynamic_cast<const DisplacementTransformType *>(compositeTransform->GetBackTransform());

      // average transformed images (starts with an all-zero image)
      typename TemplateImageType::Pointer resampledImage =
        ResampleToTarget<TemplateImageType, ImageType>(image, xavg, compositeTransform);
      // WriteImage(resampledImage, "resampledImage" + std::to_string(i) + "_" + std::to_string(k) + ".nrrd"); // debug
      xavgNew = ScaleAndAdd<TemplateImageType>(xavgNew, resampledImage, m_Weights[k]);
      // WriteImage(xavgNew, "xavgNew" + std::to_string(i) + "_" + std::to_string(k) + ".nrrd"); // debug

      if (k == 0 && dfTransform != nullptr) // transforms have a deformation field component
      {
        auto firstDF = dfTransform->GetDisplacementField();
        wavg = DisplacementImageType::New();
        wavg->CopyInformation(firstDF);
        wavg->SetRegions(firstDF->GetLargestPossibleRegion());
        wavg->Allocate(true); // initialize to zero
      }

      if (wavg != nullptr) // average the deformation fields
      {
        // average the deformation fields
        assert(wavg->IsSameImageGeometryAs(dfTransform->GetDisplacementField()));
        typename DisplacementImageType::PixelType weight;
        weight.Fill(m_Weights[k]);
        wavg = ScaleAndAdd(wavg, dfTransform->GetDisplacementField(), weight);
      }

      if (m_KeepTransforms || (compositeTransform->GetNumberOfTransforms() == 1 && affineList[k] != nullptr))
      {
        // if the composite transform is just an affine, keep it regardless of the setting
        m_TransformList[k] = const_cast<CompositeTransformType *>(compositeTransform);
      }
      this->UpdateProgress(0.01f + progressStep * (i * imageListSize + (k + 1)));
    } // for k in m_ImageList
    // WriteImage(xavgNew, "xavgNew" + std::to_string(i) + ".nrrd"); // debug

    typename AffineType::Pointer avgAffine = AffineType::New();
    typename AffineType::Pointer avgAffineInverse = AffineType::New();

    if (std::count(affineList.begin(), affineList.end(), nullptr) == 0) // all affines present
    {
      using WarperType = AverageAffineTransformFunction<AffineType>;
      WarperType average_func;
      average_func.verbose = false;
      average_func.useRigid = !m_UseNoRigid;
      for (unsigned k = 0; k < imageListSize; ++k)
      {
        average_func.PushBackAffineTransform(affineList[k], m_Weights[k]);
      }
      average_func.AverageMultipleAffineTransform(affineList[0]->GetCenter(), avgAffine);

      bool inverseExists = avgAffine->GetInverse(avgAffineInverse);
      assert(inverseExists);
      // WriteTransform(avgAffineInverse.GetPointer(), "avgAffineInverse" + std::to_string(i) + ".tfm"); // debug
    }

    if (wavg != nullptr) // we have displacement fields
    {
      // wavg *= -gradient_step
      typename DisplacementImageType::PixelType wscl;
      wscl.Fill(-1.0 * this->GetGradientStep());
      using MulType = MultiplyImageFilter<DisplacementImageType, DisplacementImageType, DisplacementImageType>;
      typename MulType::Pointer mul = MulType::New();
      mul->SetInPlace(true);
      mul->SetInput1(wavg);
      mul->SetConstant2(wscl);
      mul->Update();
      wavg = mul->GetOutput();

      typename DisplacementImageType::Pointer wavgA =
        ResampleToTarget<DisplacementImageType, DisplacementImageType>(wavg, xavgNew, avgAffineInverse);
      // WriteImage(wavgA, "wavgA" + std::to_string(i) + ".nrrd"); // debug

      typename DisplacementTransformType::Pointer wavgTransform = DisplacementTransformType::New();
      wavgTransform->SetDisplacementField(wavgA);

      typename CompositeTransformType::Pointer combinedTransform = CompositeTransformType::New();
      combinedTransform->AddTransform(avgAffineInverse);
      combinedTransform->AddTransform(wavgTransform);
      xavg = ResampleToTarget<TemplateImageType, TemplateImageType>(xavgNew, xavgNew, combinedTransform);
      // WriteTransform(combinedTransform.GetPointer(), "combinedTransform" + std::to_string(i) + ".tfm"); // debug
    }
    else // we only have the affine transforms
    {
      xavg = ResampleToTarget<TemplateImageType, TemplateImageType>(xavgNew, xavgNew, avgAffineInverse);
    }
    // WriteImage(xavg, "xavg" + std::to_string(i) + ".nrrd"); // debug

    if (m_BlendingWeight > 0)
    {
      using SharpenFilterType = LaplacianSharpeningImageFilter<TemplateImageType, TemplateImageType>;
      typename SharpenFilterType::Pointer sharpenFilter = SharpenFilterType::New();
      sharpenFilter->SetInput(xavg);
      sharpenFilter->Update();
      typename TemplateImageType::Pointer sharpened = sharpenFilter->GetOutput();
      // WriteImage(sharpened, "sharpened" + std::to_string(i) + ".nrrd"); // debug

      using AddImageFilterType = WeightedAddImageFilter<TemplateImageType, TemplateImageType, TemplateImageType>;
      typename AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
      addImageFilter->SetInPlace(true);
      addImageFilter->SetInput1(xavg);
      addImageFilter->SetInput2(sharpened);
      addImageFilter->SetAlpha(m_BlendingWeight);
      addImageFilter->Update();
      xavg = addImageFilter->GetOutput();
    }

    // WriteImage(xavg, "avgTemplate" + std::to_string(i) + ".nrrd"); // debug

    outputPtr->SetPixelContainer(xavg->GetPixelContainer());
    outputPtr->Modified();
  }

  // this->GetOutput()->Modified();  // does not help
  this->UpdateProgress(1.0);
}

} // end namespace itk

#endif // itkANTSGroupwiseBuildTemplate_hxx
