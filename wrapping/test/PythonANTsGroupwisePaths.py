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

parser = argparse.ArgumentParser(description="Python groupwise registration test.")
parser.add_argument("-i", "--input-image", action="append", help="Input image(s)", required=True)
parser.add_argument("-t", "--initial-template", required=False)
parser.add_argument("-o", "--output-template", required=True)
args = parser.parse_args()

itk.auto_progress(2)

first_image = itk.imread(args.input_image[0])

PixelType = itk.template(first_image)[1][0]
Dimension = first_image.GetImageDimension()
ImageType = itk.Image[PixelType, Dimension]
assert ImageType == type(first_image)
TemplateImageType = itk.Image[itk.F, Dimension]

gwr = itk.ANTSGroupwiseBuildTemplate[ImageType, TemplateImageType, itk.F].New()
# pairwise registration has different order of template parameters
pwr = itk.ANTSRegistration[TemplateImageType, ImageType, itk.F].New()

if args.initial_template is not None:
    print(f"Reading initial template: {args.initial_template}")
    initial_template = itk.imread(args.initial_template)
    gwr.SetInitialTemplateImage(initial_template)

gwr.SetPathList(args.input_image)
gwr.SetIterations(4)
gwr.SetGradientStep(0.2)

pwr.SetTypeOfTransform("SyN")
pwr.SetSynMetric("CC")
pwr.SetSynIterations([100, 70, 50, 0])
pwr.SetRandomSeed(30101983)  # improve reproducibility

gwr.SetPairwiseRegistration(pwr)
print(gwr)

print("Performing groupwise registration. This will take a while...")
gwr.Update()

print(f"Writing resulting template image into file: {args.output_template}")
itk.imwrite(gwr.GetOutput(), args.output_template, compression=False)

print("Test finished")
