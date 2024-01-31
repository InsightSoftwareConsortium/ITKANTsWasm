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
#ifndef itkANTSRegistration_hxx
#define itkANTSRegistration_hxx

#include "itkANTSRegistration.h"

#include <sstream>

#include "itkCastImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkPrintHelper.h"

namespace itk
{

template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::ANTSRegistration()
{
  ProcessObject::SetNumberOfRequiredOutputs(2);
  ProcessObject::SetNumberOfRequiredInputs(2);
  ProcessObject::SetNumberOfIndexedInputs(3);
  ProcessObject::SetNumberOfIndexedOutputs(2);

  SetPrimaryInputName("FixedImage");
  AddRequiredInputName("MovingImage", 1);
  AddOptionalInputName("InitialTransform", 2);
  SetPrimaryOutputName("ForwardTransform");
  // AddRequiredOutputName("InverseTransform", 1); // this method does not exist in ProcessObject

  this->ProcessObject::SetNthOutput(0, MakeOutput(0));
  this->ProcessObject::SetNthOutput(1, MakeOutput(1));
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::PrintSelf(std::ostream & os, Indent indent) const
{
  using namespace print_helper;
  Superclass::PrintSelf(os, indent);
  os << indent << "TypeOfTransform: " << this->m_TypeOfTransform << std::endl;
  os << indent << "AffineMetric: " << this->m_AffineMetric << std::endl;
  os << indent << "SynMetric: " << this->m_SynMetric << std::endl;

  os << indent << "GradientStep: " << this->m_GradientStep << std::endl;
  os << indent << "FlowSigma: " << this->m_FlowSigma << std::endl;
  os << indent << "m_TotalSigma: " << this->m_TotalSigma << std::endl;
  os << indent << "m_SamplingRate: " << this->m_SamplingRate << std::endl;
  os << indent << "NumberOfBins: " << this->m_NumberOfBins << std::endl;
  os << indent << "RandomSeed: " << this->m_RandomSeed << std::endl;
  os << indent << "SmoothingInPhysicalUnits: " << (this->m_SmoothingInPhysicalUnits ? "On" : "Off") << std::endl;

  os << indent << "SynIterations: " << this->m_SynIterations << std::endl;
  os << indent << "AffineIterations: " << this->m_AffineIterations << std::endl;
  os << indent << "ShrinkFactors: " << this->m_ShrinkFactors << std::endl;
  os << indent << "SmoothingSigmas: " << this->m_SmoothingSigmas << std::endl;

  this->m_Helper->Print(os, indent);
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::SetFixedImage(const FixedImageType * image)
{
  if (image != this->GetFixedImage())
  {
    this->ProcessObject::SetNthInput(0, const_cast<FixedImageType *>(image));
    this->Modified();
  }
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::GetFixedImage() const -> const FixedImageType *
{
  return static_cast<const FixedImageType *>(this->ProcessObject::GetInput(0));
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::SetMovingImage(const MovingImageType * image)
{
  if (image != this->GetMovingImage())
  {
    this->ProcessObject::SetNthInput(1, const_cast<MovingImageType *>(image));
    this->Modified();
  }
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::GetMovingImage() const -> const MovingImageType *
{
  return static_cast<const MovingImageType *>(this->ProcessObject::GetInput(1));
}

template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::GetWarpedMovingImage() const ->
  typename MovingImageType::Pointer
{
  using ResampleFilterType = ResampleImageFilter<MovingImageType, MovingImageType>;
  typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
  resampleFilter->SetInput(this->GetMovingImage());
  resampleFilter->SetTransform(this->GetForwardTransform());
  resampleFilter->SetSize(this->GetFixedImage()->GetLargestPossibleRegion().GetSize());
  resampleFilter->SetOutputOrigin(this->GetFixedImage()->GetOrigin());
  resampleFilter->SetOutputSpacing(this->GetFixedImage()->GetSpacing());
  resampleFilter->SetOutputDirection(this->GetFixedImage()->GetDirection());
  resampleFilter->Update();
  return resampleFilter->GetOutput();
}

template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::GetWarpedFixedImage() const ->
  typename FixedImageType::Pointer
{
  using ResampleFilterType = ResampleImageFilter<FixedImageType, FixedImageType>;
  typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
  resampleFilter->SetInput(this->GetFixedImage());
  resampleFilter->SetTransform(this->GetInverseTransform());
  resampleFilter->SetSize(this->GetMovingImage()->GetLargestPossibleRegion().GetSize());
  resampleFilter->SetOutputOrigin(this->GetMovingImage()->GetOrigin());
  resampleFilter->SetOutputSpacing(this->GetMovingImage()->GetSpacing());
  resampleFilter->SetOutputDirection(this->GetMovingImage()->GetDirection());
  resampleFilter->Update();
  return resampleFilter->GetOutput();
}

template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::SetFixedMask(const LabelImageType * mask)
{
  if (mask != this->GetFixedMask())
  {
    this->ProcessObject::SetInput("FixedMask", const_cast<LabelImageType *>(mask));
    this->Modified();
  }
}

template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::GetFixedMask() const -> const LabelImageType *
{
  return static_cast<const LabelImageType *>(this->ProcessObject::GetInput("FixedMask"));
}

template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
inline void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::SetMovingMask(const LabelImageType * mask)
{
  if (mask != this->GetMovingMask())
  {
    this->ProcessObject::SetInput("MovingMask", const_cast<LabelImageType *>(mask));
    this->Modified();
  }
}

template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::GetMovingMask() const -> const LabelImageType *
{
  return static_cast<const LabelImageType *>(this->ProcessObject::GetInput("MovingMask"));
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::GetOutput(DataObjectPointerArraySizeType index)
  -> DecoratedOutputTransformType *
{
  return static_cast<DecoratedOutputTransformType *>(this->ProcessObject::GetOutput(index));
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::GetOutput(DataObjectPointerArraySizeType index) const
  -> const DecoratedOutputTransformType *
{
  return static_cast<const DecoratedOutputTransformType *>(this->ProcessObject::GetOutput(index));
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::AllocateOutputs()
{
  const DecoratedOutputTransformType * decoratedOutputForwardTransform = this->GetOutput(0);
  if (!decoratedOutputForwardTransform || !decoratedOutputForwardTransform->Get())
  {
    this->ProcessObject::SetNthOutput(0, MakeOutput(0));
  }

  const DecoratedOutputTransformType * decoratedOutputInverseTransform = this->GetOutput(1);
  if (!decoratedOutputInverseTransform || !decoratedOutputInverseTransform->Get())
  {
    this->ProcessObject::SetNthOutput(1, MakeOutput(1));
  }
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::MakeOutput(DataObjectPointerArraySizeType)
  -> DataObjectPointer
{
  typename OutputTransformType::Pointer ptr;
  Self::MakeOutputTransform(ptr);
  typename DecoratedOutputTransformType::Pointer decoratedOutputTransform = DecoratedOutputTransformType::New();
  decoratedOutputTransform->Set(ptr);
  return decoratedOutputTransform;
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
template <typename TImage>
auto
itk::ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::CastImageToInternalType(
  const TImage * inputImage) -> typename InternalImageType::Pointer
{
  using CastFilterType = CastImageFilter<TImage, InternalImageType>;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(inputImage);
  castFilter->Update();
  typename InternalImageType::Pointer outputImage = castFilter->GetOutput();
  outputImage->DisconnectPipeline();
  return outputImage;
}

template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::GenerateData()
{
  this->AllocateOutputs();

  this->UpdateProgress(0.01);
  std::stringstream ss;
  m_Helper->SetLogStream(ss);

  const DecoratedInitialTransformType * decoratedInitialTransform = this->GetInitialTransformInput();
  if (decoratedInitialTransform != nullptr)
  {
    const InitialTransformType * initialTransform = decoratedInitialTransform->Get();
    if (initialTransform != nullptr)
    {
      m_Helper->SetMovingInitialTransform(initialTransform);
    }
  }

  typename LabelImageType::Pointer fixedMask(const_cast<LabelImageType *>(this->GetFixedMask()));
  if (fixedMask != nullptr)
  {
    m_Helper->AddFixedImageMask(fixedMask);
  }
  typename LabelImageType::Pointer movingMask(const_cast<LabelImageType *>(this->GetMovingMask()));
  if (movingMask != nullptr)
  {
    m_Helper->AddMovingImageMask(movingMask);
  }

  std::string whichTransform = this->GetTypeOfTransform();
  std::transform(whichTransform.begin(), whichTransform.end(), whichTransform.begin(), tolower);
  typename RegistrationHelperType::XfrmMethod xfrmMethod = m_Helper->StringToXfrmMethod(whichTransform);

  switch (xfrmMethod)
  {
    case RegistrationHelperType::Affine: {
      m_Helper->AddAffineTransform(m_GradientStep);
    }
    break;
    case RegistrationHelperType::Rigid: {
      m_Helper->AddRigidTransform(m_GradientStep);
    }
    break;
    case RegistrationHelperType::CompositeAffine: {
      m_Helper->AddCompositeAffineTransform(m_GradientStep);
    }
    break;
    case RegistrationHelperType::Similarity: {
      m_Helper->AddSimilarityTransform(m_GradientStep);
    }
    break;
    case RegistrationHelperType::Translation: {
      m_Helper->AddTranslationTransform(m_GradientStep);
    }
    break;
    default:
      itkExceptionMacro(<< "Unsupported transform type: " << this->GetTypeOfTransform());
  }

  // set the vector-vector parameters
  m_Helper->SetIterations({ m_AffineIterations });
  m_Helper->SetSmoothingSigmas({ m_SmoothingSigmas });
  m_Helper->SetSmoothingSigmasAreInPhysicalUnits({ m_SmoothingInPhysicalUnits });
  m_Helper->SetShrinkFactors({ m_ShrinkFactors });
  if (m_RandomSeed != 0)
  {
    m_Helper->SetRegistrationRandomSeed(m_RandomSeed);
  }

  // match the length of the iterations vector by these defaulted parameters
  std::vector<double> weights(m_AffineIterations.size(), 1.0);
  m_Helper->SetRestrictDeformationOptimizerWeights({ weights });
  std::vector<unsigned int> windows(m_AffineIterations.size(), 10);
  m_Helper->SetConvergenceWindowSizes({ windows });
  std::vector<double> thresholds(m_AffineIterations.size(), 1e-6);
  m_Helper->SetConvergenceThresholds({ thresholds });

  std::string metricType = this->GetAffineMetric();
  std::transform(metricType.begin(), metricType.end(), metricType.begin(), tolower);
  if (metricType == "jhmi")
  {
    metricType = "mi2"; // ANTs uses "mi" for Mattes MI, see:
    // https://github.com/ANTsX/ANTs/blob/v2.5.1/Examples/itkantsRegistrationHelper.hxx#L145-L152
  }
  typename RegistrationHelperType::MetricEnumeration currentMetric = m_Helper->StringToMetricType(metricType);

  typename InternalImageType::Pointer fixedImage = this->CastImageToInternalType(this->GetFixedImage());
  typename InternalImageType::Pointer movigImage = this->CastImageToInternalType(this->GetMovingImage());
  this->UpdateProgress(0.1);

  m_Helper->AddMetric(currentMetric,
                      fixedImage,
                      movigImage,
                      nullptr,
                      nullptr,
                      nullptr,
                      nullptr,
                      0u,
                      1.0,
                      RegistrationHelperType::regular,
                      m_NumberOfBins,
                      m_Radius,
                      m_UseGradientFilter,
                      false,
                      1.0,
                      50u,
                      1.1,
                      false,
                      m_SamplingRate,
                      std::sqrt(5),
                      std::sqrt(5));

  int retVal = m_Helper->DoRegistration();
  this->UpdateProgress(0.95);
  if (retVal != EXIT_SUCCESS)
  {
    itkExceptionMacro(<< "Registration failed. Helper's accumulated output:\n " << ss.str());
  }

  typename OutputTransformType::Pointer forwardTransform = m_Helper->GetModifiableCompositeTransform();
  this->SetForwardTransform(forwardTransform);

  typename OutputTransformType::Pointer inverseTransform = OutputTransformType::New();
  if (forwardTransform->GetInverse(inverseTransform))
  {
    this->SetInverseTransform(inverseTransform);
  }

  this->UpdateProgress(1.0);
}

} // end namespace itk

#endif // itkANTSRegistration_hxx
