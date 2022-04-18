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
#ifndef itkANTSRegistration_h
#define itkANTSRegistration_h

#include "itkProcessObject.h"
#include "itkImage.h"
#include "itkCompositeTransform.h"
#include "itkDataObjectDecorator.h"

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
  using MovingImageType = TTMovingImage;
  using FixedPixelType = typename InputImageType::PixelType;
  using MovingPixelType = typename OutputImageType::PixelType;

  using ParametersValueType = TParametersValueType;
  using TransformType = Transform<TParametersValueType, ImageDimension, ImageDimension>;
  using InitialTransformType = TransformType;
  using CompositeTransformType = CompositeTransform<ParametersValueType, ImageDimension>;
  using OutputTransformType = CompositeTransformType;

  /** Standard class aliases. */
  using Self = ANTSRegistration<FixedImageType, MovingImageType>;
  using Superclass = ProcessObject
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information. */
  itkTypeMacro(ANTSRegistration, ProcessObject);

  /** Standard New macro. */
  itkNewMacro(Self);

  /** Set/get the fixed image. */
  virtual void SetFixedImage(const FixedImageType * image);
  virtual const FixedImageType * GetFixedImage() const;

  /** Set/get the fixed image. */
  virtual void SetMovingImage(const MovingImageType * image);
  virtual const MovingImageType * GetMovingImage() const;

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
   * - "SyN": Symmetric normalization: Affine + deformable transformation, with mutual information as optimization metric.
   * - "SyNCC": SyN, but with cross-correlation as the metric.
   */
  itkSetStringMacro(TypeOfTransform);
  itkGetStringMacro(TypeOfTransform);

  virtual void SetInitialTransform(const InitialTransformType * transform);
  virtual const InitialTransformType * GetInitialTransform() const;

  /** Returns the transform resulting from the registration process  */
  virtual OutputTransformType *
  GetModifiableForwardTransform();
  virtual const OutputTransformType *
  GetForwardTransform() const;

  /** Returns the inverse transform resulting from the registration process, if available  */
  virtual OutputTransformType *
  GetModifiableInverseTransform();
  virtual const OutputTransformType *
  GetInverseTransform() const;

protected:
  ANTSRegistration();
  ~ANTSRegistration() override = default;

  /** Make a DataObject of the correct type to be used as the specified output. */
  using DataObjectPointerArraySizeType = ProcessObject::DataObjectPointerArraySizeType;
  using Superclass::MakeOutput;
  DataObjectPointer MakeOutput(DataObjectPointerArraySizeType) override;

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

  static void
  MakeOutputTransform(SmartPointer<InitialTransformType> & ptr)
  {
    ptr = IdentityTransform<RealType, ImageDimension>::New().GetPointer();
  }

  std::string m_TypeOfTransform{"Affine"};

private:
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
