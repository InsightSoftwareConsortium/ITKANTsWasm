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
#ifndef itkANTSRegistration_h
#define itkANTSRegistration_h

#include "itkProcessObject.h"
#include "itkImage.h"
#include "itkCompositeTransform.h"
#include "itkDataObjectDecorator.h"
#include "itkantsRegistrationHelper.h"

namespace itk
{

/** \class ANTSRegistration
 *
 * \brief Image-to-image registration method parameterized according to ANTsR or ANTsPy.
 *
 * \ingroup ANTsWasm
 *
 */
template <typename TFixedImage, typename TMovingImage, typename TParametersValueType = double>
class ANTSRegistration : public ProcessObject
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(ANTSRegistration);

  static constexpr unsigned int ImageDimension = TFixedImage::ImageDimension;

  using FixedImageType = TFixedImage;
  using MovingImageType = TMovingImage;
  using FixedPixelType = typename FixedImageType::PixelType;
  using MovingPixelType = typename MovingImageType::PixelType;

  using ParametersValueType = TParametersValueType;
  using TransformType = Transform<TParametersValueType, ImageDimension, ImageDimension>;
  using InitialTransformType = TransformType;
  using CompositeTransformType = CompositeTransform<ParametersValueType, ImageDimension>;
  using OutputTransformType = CompositeTransformType;
  using DecoratedInitialTransformType = DataObjectDecorator<InitialTransformType>;
  using DecoratedOutputTransformType = DataObjectDecorator<OutputTransformType>;

  /** Standard class aliases. */
  using Self = ANTSRegistration<FixedImageType, MovingImageType>;
  using Superclass = ProcessObject;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(ANTSRegistration, ProcessObject);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Set/get the fixed image. */
  virtual void
  SetFixedImage(const FixedImageType * image);
  virtual const FixedImageType *
  GetFixedImage() const;

  /** Set/get the fixed image. */
  virtual void
  SetMovingImage(const MovingImageType * image);
  virtual const MovingImageType *
  GetMovingImage() const;

  /** Get the moving image resampled onto fixed image grid.
   * Available after a call to Update(). Computationally expensive. */
  virtual typename MovingImageType::Pointer
  GetWarpedMovingImage() const;

  /** Get the fixed image resampled onto moving image grid.
   * Not available before a call to Update(). Computationally expensive.
   * This method raises an exception if inverse transform is not available. */
  virtual typename FixedImageType::Pointer
  GetWarpedFixedImage() const;

  /** Set the type of transformation to be optimized. A setting defines
   * a set of transformation parameterizations that are optimized,
   * the similarity metric used, and optimization parameters.
   *
   * Supported transformation types are:
   *
   * - "Translation": Translation transformation.
   * - "Rigid": Rigid transformation: Only rotation and translation.
   * - "Similarity": Similarity transformation: scaling, rotation and translation.
   * - "Affine": Affine transformation: Rigid + scaling.
   * - "AffineFast": Fast version of Affine.
   * - "SyN": Symmetric normalization: Affine + deformable transformation, with mutual information as optimization
   * metric.
   * - "SyNCC": SyN, but with cross-correlation as the metric.
   */
  itkSetStringMacro(TypeOfTransform);
  itkGetStringMacro(TypeOfTransform);

  /** Set/Get the initial transform. */
  itkSetGetDecoratedObjectInputMacro(InitialTransform, InitialTransformType);

  /** Returns the transform resulting from the registration process  */
  virtual const OutputTransformType *
  GetForwardTransform() const
  {
    return this->GetOutput(0)->Get();
  }

  /** Returns the inverse transform resulting from the registration process, if available  */
  virtual const OutputTransformType *
  GetInverseTransform() const
  {
    return this->GetOutput(1)->Get();
  }

  virtual DecoratedOutputTransformType *
  GetOutput(DataObjectPointerArraySizeType i);
  virtual const DecoratedOutputTransformType *
  GetOutput(DataObjectPointerArraySizeType i) const;

protected:
  ANTSRegistration();
  ~ANTSRegistration() override = default;

  /** Make a DataObject of the correct type to be used as the specified output. */
  using DataObjectPointerArraySizeType = ProcessObject::DataObjectPointerArraySizeType;
  using Superclass::MakeOutput;
  DataObjectPointer MakeOutput(DataObjectPointerArraySizeType) override;
  using RegistrationHelperType = ::ants::RegistrationHelper<TParametersValueType, FixedImageType::ImageDimension>;
  using InternalImageType = typename RegistrationHelperType::ImageType; // float or double pixels

  template <typename TImage>
  typename InternalImageType::Pointer
  CastImageToInternalType(const TImage *);

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData() override;

  virtual void
  AllocateOutputs();

  // helper function to create the right kind of concrete transform
  template <typename TTransform>
  static void
  MakeOutputTransform(SmartPointer<TTransform> & ptr)
  {
    ptr = TTransform::New();
  }

  /** Sets the primary output to the provided forward transform. */
  virtual void
  SetForwardTransform(const OutputTransformType * forwardTransform)
  {
    return this->GetOutput(0)->Set(forwardTransform);
  }

  /** Sets the second output to the provided inverse transform. */
  virtual void
  SetInverseTransform(const OutputTransformType * inverseTransform)
  {
    return this->GetOutput(1)->Set(inverseTransform);
  }

  std::string m_TypeOfTransform{ "Affine" };

private:
  typename RegistrationHelperType::Pointer m_Helper{ RegistrationHelperType::New() };

#ifdef ITK_USE_CONCEPT_CHECKING
  // Add concept checking such as
  // itkConceptMacro( FloatingPointPixel, ( itk::Concept::IsFloatingPoint< typename InputImageType::PixelType > ) );
#endif
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkANTSRegistration.hxx"
#endif

#endif // itkANTSRegistration
