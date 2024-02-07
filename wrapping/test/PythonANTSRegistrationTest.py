#==========================================================================
#
#   Copyright NumFOCUS
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#          https://www.apache.org/licenses/LICENSE-2.0.txt
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
#==========================================================================*/

import itk
import argparse

parser = argparse.ArgumentParser(description="Register moving image to fixed image.")
parser.add_argument("-f", "--fixed-image", required=True)
parser.add_argument("-m", "--moving-image", required=True)
parser.add_argument("-i", "--initial-transform", required=False)
parser.add_argument("-o", "--output-transform", required=True)
parser.add_argument("-r", "--resampled-moving", required=False)
args = parser.parse_args()

itk.auto_progress(2)

fixed_image = itk.imread(args.fixed_image)
moving_image = itk.imread(args.moving_image)

PixelType = itk.template(fixed_image)[1][0]
Dimension = fixed_image.GetImageDimension()
ImageType = itk.Image[PixelType, Dimension]
assert ImageType == type(fixed_image)

# registration_method = itk.ANTSRegistration[ImageType, ImageType].New(fixed_image, moving_image)
registration_method = itk.ANTSRegistration[ImageType, ImageType].New()
registration_method.SetFixedImage(fixed_image)
registration_method.SetMovingImage(moving_image)
if args.initial_transform is not None:
    transforms = itk.transformread(args.initial_transform)
    initial_transform = transforms[0]
    if len(transforms) > 1:
        print("Warning: the first transform will be used, the rest will be ignored.")
    registration_method.SetInitialTransform(initial_transform)

registration_method.Update()
forward_transform = registration_method.GetForwardTransform()
itk.transformwrite([forward_transform], args.output_transform)
inverse_transform = registration_method.GetInverseTransform()  # just check it doesn't crash

resampled_moving = itk.resample_image_filter(
    moving_image, transform=forward_transform, output_parameters_from_image=fixed_image)
itk.imwrite(resampled_moving, args.resampled_moving, compression=False)
