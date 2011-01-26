/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
#ifndef __itkModifiedGibsonLanniPointSpreadFunctionImageSource_h
#define __itkModifiedGibsonLanniPointSpreadFunctionImageSource_h

#include "itkGibsonLanniPointSpreadFunctionImageSource.h"

#include "itkGaussianImageSource.h"
#include "itkSubtractImageFilter.h"
#include "itkCommand.h"


namespace itk
{

/** \class ModifiedGibsonLanniPointSpreadFunctionImageSource
 *
 * \brief Generate a synthetic point-spread function according to the
 * Gibson-Lanni model modified by subtracting a 3D Gaussian.
 *
 * The need for this class arose from the failure of the Gibson-Lanni
 * model to adequately match a measured point-spread function.
 *
 * \ingroup DataSources Multithreaded
 */
template< class TOutputImage >
class ITK_EXPORT ModifiedGibsonLanniPointSpreadFunctionImageSource :
    public ParametricImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef ModifiedGibsonLanniPointSpreadFunctionImageSource Self;
  typedef ParametricImageSource< TOutputImage >             Superclass;
  typedef SmartPointer< Self >                              Pointer;
  typedef SmartPointer< const Self >                        ConstPointer;

  itkStaticConstMacro(ImageDimension, unsigned int,
		      TOutputImage::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ModifiedGibsonLanniPointSpreadFunctionImageSource, ParametricImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef GibsonLanniPointSpreadFunctionImageSource< TOutputImage >       GibsonLanniSourceType;
  typedef typename GibsonLanniSourceType::Pointer                         GibsonLanniSourcePointer;
  typedef GaussianImageSource< TOutputImage >                             GaussianSourceType;
  typedef typename GaussianSourceType::Pointer                            GaussianSourcePointer;
  typedef SubtractImageFilter< TOutputImage, TOutputImage, TOutputImage > SubtractImageFilterType;
  typedef typename SubtractImageFilterType::Pointer                       SubtractImageFilterPointer;

  typedef typename Superclass::ParametersValueType ParametersValueType;
  typedef typename Superclass::ParametersType      ParametersType;

  /** Set a single parameter value. */
  virtual void SetParameter(unsigned int index, ParametersValueType value);

  /** Get a single parameter value. */
  virtual ParametersValueType GetParameter(unsigned int index) const;

  /** Expects the parameters argument to contain values for ALL parameters. */
  virtual void SetParameters(const ParametersType& parameters);

  /** Gets the full parameters list. */
  virtual ParametersType GetParameters() const;

  /** Gets the total number of parameters. */
  virtual unsigned int GetNumberOfParameters() const;

  /** Callback evoked whenever the source filters are modified. */
  virtual void ModifiedCallback();

protected:
  ModifiedGibsonLanniPointSpreadFunctionImageSource();
  virtual ~ModifiedGibsonLanniPointSpreadFunctionImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();

private:
  ModifiedGibsonLanniPointSpreadFunctionImageSource(const ModifiedGibsonLanniPointSpreadFunctionImageSource&); //purposely not implemented
  void operator=(const ModifiedGibsonLanniPointSpreadFunctionImageSource&);

  GibsonLanniSourcePointer   m_GibsonLanniSource;
  GaussianSourcePointer      m_GaussianSource;
  SubtractImageFilterPointer m_SubtractFilter;

  typedef SimpleMemberCommand< Self > MemberCommandType;
  typedef typename MemberCommandType::Pointer MemberCommandPointer;
  MemberCommandPointer m_ModifiedEventCommand;
  unsigned long        m_GibsonLanniSourceObserverTag;
  unsigned long        m_GaussianSourceObserverTag;
};
} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkModifiedGibsonLanniPointSpreadFunctionImageSource.txx"
#endif

#endif
