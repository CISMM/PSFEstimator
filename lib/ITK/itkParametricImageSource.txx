/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkParametricImageSource.txx,v $
  Language:  C++
  Date:      $Date: 2009/07/17 16:10:19 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkParametricImageSource_txx
#define __itkParametricImageSource_txx
#include "itkParametricImageSource.h"

namespace itk
{

template< class TOutputImage >
ParametricImageSource< TOutputImage >
::ParametricImageSource()
{
  // Create the output. We use static_cast<> here because we know the default
  // output must be of type TOutputImage
  typename TOutputImage::Pointer output
    = static_cast<TOutputImage*>(this->MakeOutput(0).GetPointer());
  this->ProcessObject::SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput(0, output.GetPointer());

  // Set the default behavior of an image source to NOT release its
  // output bulk data prior to GenerateData() in case that bulk data
  // can be reused (an thus avoid a costly deallocate/allocate cycle).
  this->ReleaseDataBeforeUpdateFlagOff();

  // Initialize the image parameters
  this->m_Origin.Fill(0.0);
  this->m_Spacing.Fill(1.0);
  this->m_Size.Fill(1);

  // Set the transform to the identity
  this->m_Transform = TransformType::New();
}


template< class TOutputImage >
void
ParametricImageSource< TOutputImage >
::GenerateOutputInformation()
{
  OutputImageType *output;
  IndexType index = {{0}};
  SizeType size( this->m_Size );

  output = this->GetOutput(0);

  RegionType largestPossibleRegion;
  largestPossibleRegion.SetSize( size );
  largestPossibleRegion.SetIndex( index );
  output->SetLargestPossibleRegion( largestPossibleRegion );

  output->SetSpacing(this->m_Spacing);
  output->SetOrigin(this->m_Origin);
}


template< class TOutputImage >
void
ParametricImageSource< TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  unsigned int i;
  os << indent << "Origin: [";
  for ( i=0; i < this->m_Origin.Size() - 1; i++ )
    {
    os << this->m_Origin[i] << ", ";
    }
  os << this->m_Origin[i] << "]" << std::endl;

  os << indent << "Spacing: [";
  for ( i=0; i < this->m_Spacing.Size() - 1; i++ )
    {
    os << this->m_Spacing[i] << ", ";
    }
  os << this->m_Spacing[i] << "] (nanometers)" << std::endl;

  os << indent << "Size: [";
  for ( i=0; i < this->m_Size.GetSizeDimension() - 1; i++ )
    {
    os << this->m_Size[i] << ", ";
    }
  os << this->m_Size[i] << "]" << std::endl;

}

} // end namespace itk

#endif
