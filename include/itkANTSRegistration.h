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
#include "itkDisplacementFieldTransformParametersAdaptor.h"

namespace itk
{

/** \class ANTSRegistration
 *
 * \brief Image-to-image registration method parameterized according to ANTsR or ANTsPy.
 *
 * This uses image pyramids to provide reasonable registration performance and robustness.
 * Number of pyramid levels is controlled via iteration parameters.
 * There will be as many pyramid levels as there are elements in the iteration array.
 * Number of elements in shrink factor and smoothing sigma arrays (if provided) must match.
 *
 * \ingroup ANTsWasm
 * \ingroup Registration
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
  using LabelImageType = itk::Image<unsigned char, ImageDimension>;

  using ParametersValueType = TParametersValueType;
  using TransformType = Transform<TParametersValueType, ImageDimension, ImageDimension>;
  using InitialTransformType = TransformType;
  using CompositeTransformType = CompositeTransform<ParametersValueType, ImageDimension>;
  using OutputTransformType = CompositeTransformType;
  using DecoratedInitialTransformType = DataObjectDecorator<InitialTransformType>;
  using DecoratedOutputTransformType = DataObjectDecorator<OutputTransformType>;

  /** Standard class aliases. */
  using Self = ANTSRegistration<FixedImageType, MovingImageType, ParametersValueType>;
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

  /** Set/get the moving image. */
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

  /** Set/get the fixed image's mask. */
  virtual void
  SetFixedMask(const LabelImageType * mask);
  virtual const LabelImageType *
  GetFixedMask() const;

  /** Set/get the moving image's mask. */
  virtual void
  SetMovingMask(const LabelImageType * mask);
  virtual const LabelImageType *
  GetMovingMask() const;

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

  /** The metric for the affine part. Allowed metrics:
   * "MeanSquares": from MeanSquaresImageToImageMetricv4
   * "CC": neighborhood normalized cross correlation from ANTSNeighborhoodCorrelationImageToImageMetricv4
   * "GC": global normalized correlation from CorrelationImageToImageMetricv4
   * "JHMI": mutual information from JointHistogramMutualInformationImageToImageMetricv4
   * "Mattes" (default): Mattes mutual informatio from MattesMutualInformationImageToImageMetricv4
   */
  itkSetStringMacro(AffineMetric);
  itkGetStringMacro(AffineMetric);

  /** The metric for the SyN part. Allowed metrics:
   * "MeanSquares": from MeanSquaresImageToImageMetricv4
   * "CC": neighborhood normalized cross correlation from ANTSNeighborhoodCorrelationImageToImageMetricv4
   * "GC": global normalized correlation from CorrelationImageToImageMetricv4
   * "JHMI": mutual information from JointHistogramMutualInformationImageToImageMetricv4
   * "Mattes" (default): Mattes mutual informatio from MattesMutualInformationImageToImageMetricv4
   * "Demons": from DemonsImageToImageMetricv4
   */
  itkSetStringMacro(SynMetric);
  itkGetStringMacro(SynMetric);

  /** Set/Get the initial transform.
   * It transforms points from the fixed image to the moving image reference frame.
   * It is typically used to resample the moving image onto the fixed image grid. */
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

  /** Set/Get the gradient step size for transform optimizers that use it. */
  itkSetMacro(GradientStep, ParametersValueType);
  itkGetMacro(GradientStep, ParametersValueType);

  /** Set/Get smoothing for update field.
   * This only affects transform which use a deformation field. */
  itkSetMacro(FlowSigma, ParametersValueType);
  itkGetMacro(FlowSigma, ParametersValueType);

  /** Set/Get smoothing for total field.
   * This only affects transform which use a deformation field. */
  itkSetMacro(TotalSigma, ParametersValueType);
  itkGetMacro(TotalSigma, ParametersValueType);

  /** Set/Get randomg sampling percentage for estimaging the metric.
   * It is normalized to 0.0-1.0 range.
   * This can impact speed but also reproducibility and/or accuracy. */
  itkSetClampMacro(SamplingRate, ParametersValueType, 0.0, 1.0);
  itkGetMacro(SamplingRate, ParametersValueType);

  /** Set/Get number of bins for the histogram in the mutual information metric. */
  itkSetClampMacro(NumberOfBins, int, 5, NumericTraits<int>::max());
  itkGetMacro(NumberOfBins, int);

  /** Set/Get a random seed to improve reproducibility.
   * Note: ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS should be 1 for full reproducibility. */
  itkSetMacro(RandomSeed, int);
  itkGetMacro(RandomSeed, int);

  /** Set/Get radius used by neighborhood normalized cross correlation ("CC") metric. */
  itkSetMacro(Radius, unsigned int);
  itkGetMacro(Radius, unsigned int);

  /** Set/Get whether a smoothing filter is applied when estimating the metric gradient. */
  itkSetMacro(UseGradientFilter, bool);
  itkGetMacro(UseGradientFilter, bool);

  /** Set/Get whether smoothing parameters are expressed in physical units (typically millimeters). */
  itkSetMacro(SmoothingInPhysicalUnits, bool);
  itkGetMacro(SmoothingInPhysicalUnits, bool);

  /** Set/Get whether the resulting transform queue is reduced to just two transform, one linear and one deformable.
   * If false (default), there are as many transforms as there are stages. */
  itkSetMacro(CollapseCompositeTransform, bool);
  itkGetMacro(CollapseCompositeTransform, bool);

  /** Set/Get whether the all the stages should use the mask, or only the last stage (default). */
  itkSetMacro(MaskAllStages, bool);
  itkGetMacro(MaskAllStages, bool);

  /** Set/Get number of iterations for each pyramid level for SyN transforms.
   * Shrink factors and smoothing sigmas for SyN are determined based on iterations. */
  itkSetMacro(SynIterations, std::vector<unsigned int>);
  itkGetConstReferenceMacro(SynIterations, std::vector<unsigned int>);

  /** Set/Get number of iterations for each pyramid level for low-dim transforms.
   * Low dimensionality transforms are Translation, Rigid, Similarity, and Affine. */
  itkSetMacro(AffineIterations, std::vector<unsigned int>);
  itkGetConstReferenceMacro(AffineIterations, std::vector<unsigned int>);

  /** Set/Get shrink factor for each pyramid level for low-dim transforms.
   * Low dimensionality transforms are Translation, Rigid, Similarity, and Affine. */
  itkSetMacro(ShrinkFactors, std::vector<unsigned int>);
  itkGetConstReferenceMacro(ShrinkFactors, std::vector<unsigned int>);

  /** Set/Get smoothing sigmas for each pyramid level for low-dim transforms.
   * Low dimensionality transforms are Translation, Rigid, Similarity, and Affine. */
  itkSetMacro(SmoothingSigmas, std::vector<float>);
  itkGetConstReferenceMacro(SmoothingSigmas, std::vector<float>);

  /** Set/Get the optimizer weights. When set, this allows restricting the optimization
   * of the displacement field, translation, rigid or affine transform on a per-component basis.
   * For example, to limit the deformation or rotation of 3-D volume to the first two dimensions,
   * specify a weight vector of ‘(1,1,0)’ for a 3D displacement field
   * or ‘(1,1,0,1,1,0)’ for a rigid transformation. */
  itkSetMacro(RestrictTransformation, std::vector<ParametersValueType>);
  itkGetConstReferenceMacro(RestrictTransformation, std::vector<ParametersValueType>);

  /** Set/Get the subsampling factor for displacement fields results.
   * A factor of 1 results in no subsampling. This is applied in all dimensions.
   * The default is 2. */
  itkSetMacro(DisplacementFieldSubsamplingFactor, unsigned int);
  itkGetMacro(DisplacementFieldSubsamplingFactor, unsigned int);

  virtual DecoratedOutputTransformType *
  GetOutput(DataObjectPointerArraySizeType i);
  virtual const DecoratedOutputTransformType *
  GetOutput(DataObjectPointerArraySizeType i) const;
  virtual DecoratedOutputTransformType *
  GetOutput()
  {
    return this->GetOutput(0); // return the forward transform
  }
  virtual const DecoratedOutputTransformType *
  GetOutput() const
  {
    return this->GetOutput(0); // return the forward transform
  }

  using ProcessObject::AddInput;
  using ProcessObject::RemoveInput;
  using ProcessObject::GetInput;
  using ProcessObject::SetInput;

  /** Set/Get the image input of this process object.  */
  virtual void
  SetInput(const FixedImageType * input)
  {
    this->SetInput(0, input);
  }

  virtual void
  SetInput(unsigned index, const FixedImageType * image);


protected:
  ANTSRegistration();
  ~ANTSRegistration() override = default;

  /** Make a DataObject of the correct type to be used as the specified output. */
  using DataObjectPointerArraySizeType = ProcessObject::DataObjectPointerArraySizeType;
  using Superclass::MakeOutput;
  DataObjectPointer MakeOutput(DataObjectPointerArraySizeType) override;
  using RegistrationHelperType = ::ants::RegistrationHelper<TParametersValueType, FixedImageType::ImageDimension>;
  using InternalImageType = typename RegistrationHelperType::ImageType; // float or double pixels
  using DisplacementFieldTransformType = typename RegistrationHelperType::DisplacementFieldTransformType;
  using DisplacementFieldType = typename DisplacementFieldTransformType::DisplacementFieldType;
  using DisplacementFieldTransformParametersAdaptorType =
    DisplacementFieldTransformParametersAdaptor<DisplacementFieldTransformType>;

  template <typename TImage>
  typename InternalImageType::Pointer
  CastImageToInternalType(const TImage *);

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Returns true if registration was successful. */
  void
  SingleStageRegistration(typename RegistrationHelperType::XfrmMethod xfrmMethod,
                          const InitialTransformType *                initialTransform,
                          typename InternalImageType::Pointer         fixedImage,
                          typename InternalImageType::Pointer         movingImage,
                          bool                                        useMasks,
                          unsigned                                    nTimeSteps = 4);

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
  std::string m_AffineMetric{ "Mattes" };
  std::string m_SynMetric{ "Mattes" };

  ParametersValueType m_GradientStep{ 0.2 };
  ParametersValueType m_FlowSigma{ 3.0 };
  ParametersValueType m_TotalSigma{ 0.0 };
  ParametersValueType m_SamplingRate{ 0.2 };
  int                 m_NumberOfBins{ 32 };
  int                 m_RandomSeed{ 0 };
  bool                m_SmoothingInPhysicalUnits{ false };
  bool                m_UseGradientFilter{ false };
  unsigned int        m_Radius{ 4 };
  bool                m_CollapseCompositeTransform{ true };
  bool                m_MaskAllStages{ false };
  unsigned int        m_DisplacementFieldSubsamplingFactor{ 2 };

  std::vector<unsigned int> m_SynIterations{ 40, 20, 0 };
  std::vector<unsigned int> m_AffineIterations{ 2100, 1200, 1200, 10 };
  std::vector<unsigned int> m_ShrinkFactors{ 6, 4, 2, 1 };
  std::vector<float>        m_SmoothingSigmas{ 3, 2, 1, 0 };

  std::vector<ParametersValueType> m_RestrictTransformation;

private:
  typename RegistrationHelperType::Pointer                          m_Helper{ RegistrationHelperType::New() };
  typename DisplacementFieldTransformParametersAdaptorType::Pointer m_DisplacementFieldAdaptor{
    DisplacementFieldTransformParametersAdaptorType::New()
  };

#ifdef ITK_USE_CONCEPT_CHECKING
  static_assert(TFixedImage::ImageDimension == TMovingImage::ImageDimension,
                "Fixed and moving images must have the same dimension.");
  static_assert(ImageDimension >= 2, "Images must be at least two-dimensional.");
#endif
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkANTSRegistration.hxx"
#endif

#endif // itkANTSRegistration
