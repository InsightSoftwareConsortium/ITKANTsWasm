/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
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

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{

template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>
::ANTSRegistration()
{
  ProcessObject::SetNumberOfRequiredOutputs(2);
  ProcessObject::SetNumberOfRequiredInputs(2);
  ProcessObject::SetNumberOfIndexedInputs(3);
  ProcessObject::SetNumberOfIndexedOutputs(2);
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "TypeOfTransform: " << this->m_TypeOfTransform << std::endl;
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>
::SetFixedImage(const FixedImageType * image)
{
  if (image != static_cast<FixedImageType *>(this->ProcessObject::GetInput(0)))
  {
    this->ProcessObject::SetNthInput(0, const_cast<FixedImageType *>(image));
    this->Modified();
  }
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>
::GetFixedImage() const -> const FixedImageType *
{
  return static_cast<const FixedImageType *>(this->ProcessObject::GetInput(0));
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>
::SetMovingImage(const MovingImageType * image)
{
  if (image != static_cast<MovingImageType *>(this->ProcessObject::GetInput(1)))
  {
    this->ProcessObject::SetNthInput(1, const_cast<MovingImageType *>(image));
    this->Modified();
  }
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>
::GetMovingImage() const -> const MovingImageType *
{
  return static_cast<const MovingImageType *>(this->ProcessObject::GetInput(1));
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>
::GetOutput(DataObjectPointerArraySizeType index) -> DecoratedOutputTransformType *
{
  return static_cast<DecoratedOutputTransformType *>(this->ProcessObject::GetOutput(index));
}

template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
auto
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>
::GetOutput(DataObjectPointerArraySizeType index) const -> const DecoratedOutputTransformType *
{
  return static_cast<const DecoratedOutputTransformType *>(this->ProcessObject::GetOutput(index));
}

template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>
::AllocateOutputs()
{
  DecoratedOutputTransformType * decoratedOutputForwardTransform = this->GetOutputForwardTransform();
  if (!decoratedOutputForwardTransform->Get())
  {
    typename OutputTransformType::Pointer ptr;
    Self::MakeOutputTransform(ptr);
    decoratedOutputForwardTransform->Set(ptr);
  }

  DecoratedOutputTransformType * decoratedOutputInverseTransform = this->GetOutputInverseTransform();
  if (!decoratedOutputInverseTransform->Get())
  {
    typename OutputTransformType::Pointer ptr;
    Self::MakeOutputTransform(ptr);
    decoratedOutputInverseTransform->Set(ptr);
  }
}


template <typename TFixedImage, typename TMovingImage, typename TParametersValueType>
void
ANTSRegistration<TFixedImage, TMovingImage, TParametersValueType>
::GenerateData()
{
  this->AllocateOutputs();

}

} // end namespace itk

#endif // itkANTSRegistration_hxx
