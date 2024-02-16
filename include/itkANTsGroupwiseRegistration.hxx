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
#ifndef itkANTsGroupwiseRegistration_hxx
#define itkANTsGroupwiseRegistration_hxx

#include <sstream>

#include "itkPrintHelper.h"
#include "itkANTsGroupwiseRegistration.h" // needed by VS code completion
#include "itkImageDuplicator.h"
#include "itkWeightedAddImageFilter.h"

namespace itk
{

template <typename TImage, typename TTemplateImage, typename TParametersValueType>
ANTsGroupwiseRegistration<TImage, TTemplateImage, TParametersValueType>::ANTsGroupwiseRegistration()
{
  this->SetPrimaryInputName("InitialTemplate");
  this->SetPrimaryOutputName("OptimizedImage");

  this->SetInput(0, TemplateImageType::New()); // empty initial template

  this->GetMultiThreader()->SetMaximumNumberOfThreads(1); // registrations are already multi-threaded
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
void
ANTsGroupwiseRegistration<TImage, TTemplateImage, TParametersValueType>::PrintSelf(std::ostream & os,
                                                                                   Indent         indent) const
{
  using namespace print_helper;
  Superclass::PrintSelf(os, indent);

  os << indent << "GradientStep: " << this->m_GradientStep << std::endl;
  os << indent << "BlendingWeight: " << this->m_BlendingWeight << std::endl;
  os << indent << "UseNoRigid: " << (this->m_UseNoRigid ? "On" : "Off") << std::endl;
  os << indent << "Iterations: " << this->m_Iterations << std::endl;
  os << indent << "Weights: " << this->m_Weights << std::endl;

  os << indent << "ImageList: " << std::endl;
  unsigned i = 0;
  for (const auto & image : this->m_ImageList)
  {
    os << indent.GetNextIndent() << "Image" << i++ << ": ";
    image->Print(os, indent.GetNextIndent());
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
ANTsGroupwiseRegistration<TImage, TTemplateImage, TParametersValueType>::VerifyInputInformation() const
{
  if (m_ImageList.empty())
  {
    itkExceptionMacro("No input images provided.");
  }

  if (m_Weights.size() > 0 && m_Weights.size() != m_ImageList.size())
  {
    itkExceptionMacro("The number of weights is different from the number of images.");
  }

  if (m_ImageList.size() < 2)
  {
    itkExceptionMacro("At least two input images are required.");
  }
}

template <typename TImage, typename TTemplateImage, typename TParametersValueType>
DataObject::Pointer
ANTsGroupwiseRegistration<TImage, TTemplateImage, TParametersValueType>::MakeOutput(DataObjectPointerArraySizeType)
{
  typename OutputTransformType::Pointer ptr;
  Self::MakeOutputTransform(ptr);
  typename DecoratedOutputTransformType::Pointer decoratedOutputTransform = DecoratedOutputTransformType::New();
  decoratedOutputTransform->Set(ptr);
  return decoratedOutputTransform;
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
auto
ANTsGroupwiseRegistration<TImage, TTemplateImage, TParametersValueType>::DuplicateImage(const TemplateImageType * image)
  -> typename TemplateImageType::Pointer
{
  using DuplicatorType = ImageDuplicator<TemplateImageType>;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(image);
  duplicator->Update();
  return duplicator->GetOutput();
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
auto
ANTsGroupwiseRegistration<TImage, TTemplateImage, TParametersValueType>::ResampleToTarget(
  const ImageType *                    input,
  const TemplateImageType *            target,
  typename TransformType::ConstPointer transform) -> typename TemplateImageType::Pointer
{
  using ResampleFilterType =
    ResampleImageFilter<ImageType, TemplateImageType, ParametersValueType, ParametersValueType>;
  typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
  resampleFilter->SetInput(input);
  if (transform)
  {
    resampleFilter->SetTransform(transform);
  }
  resampleFilter->SetOutputParametersFromImage(target);
  resampleFilter->Update();
  return resampleFilter->GetOutput();
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
template <typename OutImageType, typename ArrayImageType>
auto
ANTsGroupwiseRegistration<TImage, TTemplateImage, TParametersValueType>::AverageTransformedImages(
  const std::vector<typename ArrayImageType::Pointer> &  imageList,
  const std::vector<typename AffineType::ConstPointer> & affineList) -> typename OutImageType::Pointer
{
  assert(imageList.size() == affineList.size());
  assert(imageList.size() == m_Weights.size());

  typename OutImageType::Pointer average = OutImageType::New();
  average->CopyInformation(imageList[0]);
  average->SetRegions(imageList[0]->GetLargestPossibleRegion());
  average->Allocate(true); // initialize to zero

  using AddImageFilterType = WeightedAddImageFilter<OutImageType, OutImageType, OutImageType>;
  typename AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
  for (unsigned i = 0; i < imageList.size(); ++i)
  {
    typename OutImageType::Pointer resampledImage;
    if (affineList[i] ||
        !average->SameImageGridAs(imageList[i], this->GetCoordinateTolerance(), this->GetDirectionTolerance()))
    {
      resampledImage = ResampleToTarget(imageList[i], average, affineList[i]);
    }
    else // no need to resample (should be the case with displacement fields)
    {
      if constexpr (std::is_same_v<OutImageType, ImageType>)
      {
        resampledImage = imageList[i];
      }
      else // we need to cast
      {
        using CastFilterType = CastImageFilter<ImageType, OutImageType>;
        typename CastFilterType::Pointer castFilter = CastFilterType::New();
        castFilter->SetInput(imageList[i]);
        castFilter->Update();
        resampledImage = castFilter->GetOutput();
      }
    }

    addImageFilter->SetInput1(average);
    addImageFilter->SetInput2(resampledImage);
    addImageFilter->SetAlpha(1.0 - m_Weights[i]);
    // addImageFilter->SetInPlace(true); // breaks the result in second iteration
    addImageFilter->Update();
    average = addImageFilter->GetOutput();
  }

  return average;
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
void
ANTsGroupwiseRegistration<TImage, TTemplateImage, TParametersValueType>::GenerateData()
{
  this->UpdateProgress(0.0);

  if (!m_PairwiseRegistration) // TODO: allow setting a custom pairwise registration
  {
    m_PairwiseRegistration = PairwiseType::New();
    m_PairwiseRegistration->SetTypeOfTransform("SyN");
    m_PairwiseRegistration->SetTypeOfTransform("QuickRigid"); // debug
  }

  if (m_Weights.empty())
  {
    m_Weights.resize(m_ImageList.size(), 1.0 / m_ImageList.size());
  }
  else // normalize to sum to 1
  {
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

  typename TemplateImageType::Pointer initialTemplate = dynamic_cast<TemplateImageType *>(this->GetInput(0));
  if (!initialTemplate)
  {
    itkExceptionMacro("Initial template must be a float-pixel image.");
  }

  std::vector<typename AffineType::ConstPointer> emptyAffines(m_ImageList.size(), nullptr);
  if (initialTemplate->GetLargestPossibleRegion().GetNumberOfPixels() > 0)
  {
    // copy the initial template into the output (the current average template)
    this->GraftOutput(this->DuplicateImage(initialTemplate));
  }
  else
  {
    // create a simple average of all the input images
    this->GraftOutput(this->AverageTransformedImages<TemplateImageType, ImageType>(m_ImageList, emptyAffines));
  }
  this->UpdateProgress(0.01);

  using DisplacementTransformType = DisplacementFieldTransform<ParametersValueType, ImageDimension>;
  using DisplacementImageType = typename DisplacementTransformType::DisplacementFieldType;

  for (unsigned i = 0; i < m_Iterations; ++i)
  {
    typename TemplateImageType::Pointer                  currentTemplate = this->DuplicateImage(this->GetOutput());
    std::vector<typename AffineType::ConstPointer>       affineList(m_ImageList.size(), nullptr);
    std::vector<typename DisplacementImageType::Pointer> dfList(m_ImageList.size(), nullptr);
    for (unsigned k = 0; k < m_ImageList.size(); ++k)
    {
      m_PairwiseRegistration->SetFixedImage(currentTemplate);
      m_PairwiseRegistration->SetMovingImage(m_ImageList[k]);
      m_PairwiseRegistration->Update();

      const CompositeTransformType * compositeTransform = m_PairwiseRegistration->GetForwardTransform();

      affineList[k] = dynamic_cast<const AffineType *>(compositeTransform->GetBackTransform());
      auto dfTransform = dynamic_cast<const DisplacementTransformType *>(compositeTransform->GetFrontTransform());
      auto dfNonConstTransform = const_cast<DisplacementTransformType *>(dfTransform);
      dfList[k] = dfNonConstTransform->GetDisplacementField();

      // average transformed images, start with an all-zero image
      typename TemplateImageType::Pointer xavgNew =
        this->AverageTransformedImages<TemplateImageType, ImageType>(m_ImageList, affineList);

      // average the deformation fields
      typename DisplacementImageType::Pointer wavg =
        this->AverageTransformedImages<DisplacementImageType, DisplacementImageType>(dfList, emptyAffines);
    }

    if (m_UseNoRigid)
    {
      // avgaffine = utils.average_affine_transform_no_rigid(affineList)
    }
    else
    {
      // avgaffine = utils.average_affine_transform(affineList)
    }

    this->GraftOutput(currentTemplate);
  }

  this->UpdateProgress(0.99);

  // typename OutputTransformType::Pointer inverseTransform = OutputTransformType::New();
  // if (forwardTransform->GetInverse(inverseTransform))
  // {
  //   this->SetInverseTransform(inverseTransform);
  // }
  // else
  // {
  //   this->SetInverseTransform(nullptr);
  // }

  this->UpdateProgress(1.0);
}

} // end namespace itk

#endif // itkANTsGroupwiseRegistration_hxx
