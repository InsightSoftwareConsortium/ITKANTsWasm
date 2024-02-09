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
ANTsGroupwiseRegistration<TImage, TTemplateImage, TParametersValueType>::PrintSelf(std::ostream & os, Indent indent) const
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
auto
ANTsGroupwiseRegistration<TImage, TTemplateImage, TParametersValueType>::MakeOutput(DataObjectPointerArraySizeType) -> DataObjectPointer
{
  typename OutputTransformType::Pointer ptr;
  Self::MakeOutputTransform(ptr);
  typename DecoratedOutputTransformType::Pointer decoratedOutputTransform = DecoratedOutputTransformType::New();
  decoratedOutputTransform->Set(ptr);
  return decoratedOutputTransform;
}


template <typename TImage, typename TTemplateImage, typename TParametersValueType>
void
ANTsGroupwiseRegistration<TImage, TTemplateImage, TParametersValueType>::GenerateData()
{
  this->AllocateOutputs();

  this->UpdateProgress(0.01);

  // TODO: reimplement stuff from:
  // https://github.com/ANTsX/ANTsPy/blob/master/ants/registration/build_template.py

  this->UpdateProgress(0.95);


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
