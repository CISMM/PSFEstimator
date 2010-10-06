/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkRotationalExtrusionTransform.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkRotationalExtrusionTransform_h
#define __itkRotationalExtrusionTransform_h

#include <iostream>

#include "itkTransform.h"
#include "itkExceptionObject.h"
#include "itkMacro.h"

namespace itk
{
/**
 * RotationalExtrusion transformation used to sweep a plane from a 3D
 * image about an axis to form a 3D image.
 *
 * ScalarT       The type to be used for scalar numeric values.  Either
 *               float or double.
 *
 * \ingroup Transforms
 *
 */

template< class TScalarType = double >
// Number of dimensions in the input space
class RotationalExtrusionTransform:
    public Transform< TScalarType, 3, 3 >
{
public:
  /** Standard typedefs   */
  typedef RotationalExtrusionTransform                        Self;
  typedef Transform< TScalarType, 3, 3 >  Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(RotationalExtrusionTransform, Transform);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

  /** Dimension of the domain space. */
  itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(SpaceDimension, unsigned int, 3);
  itkStaticConstMacro( ParametersDimension, unsigned int, 0);

  /** Parameters Type   */
  typedef typename Superclass::ParametersType            ParametersType;
  typedef typename Superclass::JacobianType              JacobianType;
  typedef typename Superclass::ScalarType                ScalarType;
  typedef typename Superclass::InputPointType            InputPointType;
  typedef typename Superclass::OutputPointType           OutputPointType;
  typedef typename Superclass::InputVectorType           InputVectorType;
  typedef typename Superclass::OutputVectorType          OutputVectorType;
  typedef typename Superclass::InputVnlVectorType        InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType       OutputVnlVectorType;
  typedef typename Superclass::InputCovariantVectorType  InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType OutputCovariantVectorType;

  /** Base inverse transform type. This type should not be changed to the
   * concrete inverse transform type or inheritance would be lost.*/
  typedef typename Superclass::InverseTransformBaseType InverseTransformBaseType;
  typedef typename InverseTransformBaseType::Pointer    InverseTransformBasePointer;

  virtual OutputPointType TransformPoint(const InputPointType & inPt ) const
  {
    // Get distance from axis (in this case, the z axis at x=0, y=0
    TScalarType distance = sqrt(inPt[0]*inPt[0] + inPt[1]*inPt[1]);

    OutputPointType outPt;
    outPt[0] = distance;
    outPt[1] = 0.0;
    outPt[2] = inPt[2];

    return outPt;
  }

protected:
  /** Construct a RotationalExtrusionTransform object
   *
   * This method constructs a new RotationalExtrusionTransform object.
   */
  RotationalExtrusionTransform() : Superclass(3, 0) {};

  /** Destroy an RotationalExtrusionTransform object   */
  virtual ~RotationalExtrusionTransform() {};

  /** Print contents of an RotationalExtrusionTransform */
  void PrintSelf(std::ostream & s, Indent indent) const
  {
    Superclass::PrintSelf(s, indent);
  }

private:

  RotationalExtrusionTransform(const Self & other);
  const Self & operator=(const Self &);
}; //class RotationalExtrusionTransform
} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_RotationalExtrusionTransform(_, EXPORT, TypeX, TypeY)                \
  namespace itk                                                              \
  {                                                                          \
  _( 2 ( class EXPORT RotationalExtrusionTransform< ITK_TEMPLATE_2 TypeX > ) )            \
  namespace Templates                                                        \
  {                                                                          \
  typedef RotationalExtrusionTransform< ITK_TEMPLATE_2 TypeX >  RotationalExtrusionTransform##TypeY; \
  }                                                                          \
  }

#if ITK_TEMPLATE_EXPLICIT
#include "Templates/itkRotationalExtrusionTransform+-.h"
#endif

#endif /* __itkRotationalExtrusionTransform_h */
