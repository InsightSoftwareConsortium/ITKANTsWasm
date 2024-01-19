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

namespace itk
{

template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::ANTSRegistration()
{
  ProcessObject::SetNumberOfRequiredOutputs(2);
  ProcessObject::SetNumberOfRequiredInputs(2);
  ProcessObject::SetNumberOfIndexedInputs(3);
  ProcessObject::SetNumberOfIndexedOutputs(2);

  typename OutputTransformType::Pointer ptr;
  Self::MakeOutputTransform(ptr);
  typename DecoratedOutputTransformType::Pointer decoratedOutputTransform = DecoratedOutputTransformType::New();
  decoratedOutputTransform->Set(ptr);
  this->ProcessObject::SetNthOutput(0, decoratedOutputTransform);

  // this->SetForwardTransformInput(decoratedOutputTransform);
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "TypeOfTransform: " << this->m_TypeOfTransform << std::endl;
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::SetFixedImage(const FixedImageType * image)
{
  if (image != static_cast<FixedImageType *>(this->ProcessObject::GetInput(0)))
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
  if (image != static_cast<MovingImageType *>(this->ProcessObject::GetInput(1)))
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
  const DecoratedOutputTransformType * decoratedOutputForwardTransform = this->GetForwardTransformInput();
  if (!decoratedOutputForwardTransform->Get())
  {
    typename OutputTransformType::Pointer ptr;
    Self::MakeOutputTransform(ptr);
    typename DecoratedOutputTransformType::Pointer decoratedOutputTransform = DecoratedOutputTransformType::New();
    decoratedOutputTransform->Set(ptr);
    this->SetForwardTransformInput(decoratedOutputTransform);
  }

  const DecoratedOutputTransformType * decoratedOutputInverseTransform = this->GetInverseTransformInput();
  if (!decoratedOutputInverseTransform->Get())
  {
    typename OutputTransformType::Pointer ptr;
    Self::MakeOutputTransform(ptr);
    typename DecoratedOutputTransformType::Pointer decoratedOutputTransform = DecoratedOutputTransformType::New();
    decoratedOutputTransform->Set(ptr);
    this->SetInverseTransformInput(decoratedOutputTransform);
  }
}


// template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
// auto
// ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>::MakeOutput(DataObjectPointerArraySizeType)
//   -> DataObjectPointer
//{
//   return DataObjectPointer();
// }


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
  // this->AllocateOutputs();

  std::stringstream ss;
  m_Helper->SetLogStream(ss);

  std::string whichTransform = this->GetTypeOfTransform();
  std::transform(whichTransform.begin(), whichTransform.end(), whichTransform.begin(), tolower);
  typename RegistrationHelperType::XfrmMethod xfrmMethod = m_Helper->StringToXfrmMethod(whichTransform);

  double learningRate = 0.2; // TODO: Make this a parameter

  switch (xfrmMethod)
  {
    case RegistrationHelperType::Affine: {
      m_Helper->AddAffineTransform(learningRate);
    }
    break;
    case RegistrationHelperType::Rigid: {
      m_Helper->AddRigidTransform(learningRate);
    }
    break;
    case RegistrationHelperType::CompositeAffine: {
      m_Helper->AddCompositeAffineTransform(learningRate);
    }
    break;
    case RegistrationHelperType::Similarity: {
      m_Helper->AddSimilarityTransform(learningRate);
    }
    break;
    case RegistrationHelperType::Translation: {
      m_Helper->AddTranslationTransform(learningRate);
    }
    break;
    default:
      itkExceptionMacro(<< "Unsupported transform type: " << whichTransform);
  }

  typename RegistrationHelperType::MetricEnumeration currentMetric = m_Helper->StringToMetricType("mi");

  // assign default image metric variables
  typename RegistrationHelperType::SamplingStrategy samplingStrategy = RegistrationHelperType::none;
  int                                               numberOfBins = 32;
  unsigned int                                      radius = 4;
  bool                                              useGradientFilter = false;

  typename InternalImageType::Pointer fixedImage = this->CastImageToInternalType(this->GetFixedImage());
  typename InternalImageType::Pointer movigImage = this->CastImageToInternalType(this->GetMovingImage());

  m_Helper->AddMetric(currentMetric,
                      fixedImage,
                      movigImage,
                      nullptr,
                      nullptr,
                      nullptr,
                      nullptr,
                      0u,
                      1.0,
                      samplingStrategy,
                      numberOfBins,
                      radius,
                      useGradientFilter,
                      false,
                      1.0,
                      50u,
                      1.1,
                      false,
                      0.1,
                      std::sqrt(5),
                      std::sqrt(5));

  int retVal = m_Helper->DoRegistration();
  if (retVal != EXIT_SUCCESS)
  {
    itkExceptionMacro(<< "Registration failed. Helper's accumulated output:\n " << ss.str());
  }

  typename DecoratedOutputTransformType::Pointer decoratedOutputTransform = DecoratedOutputTransformType::New();
  decoratedOutputTransform->Set(m_Helper->GetModifiableCompositeTransform());
  this->SetForwardTransformInput(decoratedOutputTransform);
}

} // end namespace itk

#endif // itkANTSRegistration_hxx
