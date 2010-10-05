/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScanImageFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2009/09/08 21:33:37 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkScanImageFilter_cxx
#define __itkScanImageFilter_cxx

#include <itkImageLinearConstIteratorWithIndex.h>
#include <itkProgressReporter.h>
#include <itkScanImageFilter.h>

namespace itk {

/**
 * Constructor.
 */
template <class TInputImage, class TOutputImage, class TAccumulator>
ScanImageFilter<TInputImage,TOutputImage,TAccumulator>
::ScanImageFilter()
{
  this->SetNumberOfRequiredInputs(1);
  m_ScanDimension = InputImageDimension-1;
  m_ScanOrder     = INCREASING_ORDER;
}


//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage, class TAccumulator>
void
ScanImageFilter<TInputImage,TOutputImage,TAccumulator>
::GenerateInputRequestedRegion()
{
  if ( this->GetInput() )
    {
    typename TInputImage::RegionType inputRequestedRegion;
    this->GenerateInputRequestedRegionForOutputRequestedRegion
      (this->GetOutput()->GetRequestedRegion(), inputRequestedRegion);
    InputImagePointer input = const_cast< TInputImage * >( this->GetInput() );
    input->SetRequestedRegion(inputRequestedRegion);
    }
}


//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage, class TAccumulator>
void
ScanImageFilter<TInputImage,TOutputImage,TAccumulator>
::GenerateInputRequestedRegionForOutputRequestedRegion
(const OutputImageRegionType& outputRegion,
 InputImageRegionType& inputRegion)
{

  if (this->GetInput())
    {
    //typename TInputImage::RegionType requestedRegion;
    typename TInputImage::IndexType  inputIndex;
    typename TInputImage::SizeType   inputSize;
    typename TInputImage::IndexType  inputLargestIndex;
    typename TInputImage::SizeType   inputLargestSize;
    typename TOutputImage::IndexType outputIndex;
    typename TOutputImage::SizeType  outputSize;

    outputSize = outputRegion.GetSize();
    outputIndex = outputRegion.GetIndex();
    inputLargestSize  = this->GetInput()->GetLargestPossibleRegion().GetSize();
    inputLargestIndex = this->GetInput()->GetLargestPossibleRegion().GetIndex();

    for (unsigned int i = 0; i < InputImageDimension; i++)
      {
      if (i == m_ScanDimension)
        {
	if (m_ScanOrder == INCREASING_ORDER)
          {
	  inputIndex[i] = 0;
	  inputSize[i]  = outputSize[i] + outputIndex[i];
          }
        else
          {
          // DECREASING_ORDER
	  inputIndex[i] = outputIndex[i];
	  inputSize[i]  = inputLargestSize[i] - outputIndex[i];
          }
        }
      else
        {
	inputIndex[i] = outputIndex[i];
	inputSize[i]  = outputSize[i];
        }
      inputIndex[i] = outputIndex[i];
      }
    inputRegion.SetSize(inputSize);
    inputRegion.SetIndex(inputIndex);
    }
}


//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage, class TAccumulator>
int
ScanImageFilter<TInputImage,TOutputImage,TAccumulator>
::SplitRequestedRegion(int i, int num, OutputImageRegionType& splitRegion)
{
  // Get the output pointer
  OutputImageType * outputPtr = this->GetOutput();
  const typename TOutputImage::SizeType& requestedRegionSize
    = outputPtr->GetRequestedRegion().GetSize();

  int splitAxis;
  typename TOutputImage::IndexType splitIndex;
  typename TOutputImage::SizeType splitSize;

  // Initialize the splitRegion to the output requested region
  splitRegion = outputPtr->GetRequestedRegion();
  splitIndex  = splitRegion.GetIndex();
  splitSize   = splitRegion.GetSize();

  // Split on the outermost dimension available that is not equal to the
  // scan dimension.
  for (unsigned int j = outputPtr->GetImageDimension()-1; j >= 0; j--)
    {
    if (j != m_ScanDimension)
      {
      splitAxis = j;
      break;
      }
    }
  while (requestedRegionSize[splitAxis] == 1)
    {
    --splitAxis;
    if (splitAxis < 0)
      { // cannot split
      itkDebugMacro("  Cannot Split");
      return 1;
      }
    }

  // determine the actual number of pieces that will be generated
  typename TOutputImage::SizeType::SizeValueType range
    = requestedRegionSize[splitAxis];
  int valuesPerThread = (int)::vcl_ceil(range/(double)num);
  int maxThreadIdUsed = (int)::vcl_ceil(range/(double)valuesPerThread) - 1;

  // Split the region
  if (i < maxThreadIdUsed)
    {
    splitIndex[splitAxis] += i*valuesPerThread;
    splitSize[splitAxis] = valuesPerThread;
    }
  if (i == maxThreadIdUsed)
    {
    splitIndex[splitAxis] += i*valuesPerThread;
    // last thread needs to process the "rest" dimension being split
    splitSize[splitAxis] = splitSize[splitAxis] - i*valuesPerThread;
    }

  // set the split region ivars
  splitRegion.SetIndex( splitIndex );
  splitRegion.SetSize( splitSize );

  itkDebugMacro("  Split Piece: " << splitRegion );

  return maxThreadIdUsed + 1;
}


//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage, class TAccumulator>
void
ScanImageFilter<TInputImage,TOutputImage,TAccumulator>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
		       int threadId)
{
  // Use the output image to report the progress. This should be set
  // to the number of lines.
  ProgressReporter progress(this, threadId,
                            outputRegionForThread.GetNumberOfPixels());

  // Get some values, just to be easier to manipulate.
  InputImageConstPointer inputImage = this->GetInput();
  OutputImagePointer    outputImage = this->GetOutput();

  InputImageRegionType inputRegionForThread;
  this->GenerateInputRequestedRegionForOutputRequestedRegion
    (outputRegionForThread, inputRegionForThread);

  // Create the iterators for input and output image.
  typedef ImageLinearConstIteratorWithIndex<InputImageType> InputIteratorType;
  InputIteratorType iIt( inputImage, inputRegionForThread );

  iIt.SetDirection( m_ScanDimension );

  // Instantiate the accumulator.
  AccumulatorType accumulator = this->NewAccumulator(1);

  iIt.GoToBegin();

  // OK, everything is ready... lets the linear iterator do its job !
  while( !iIt.IsAtEnd() )
    {

    // Init the accumulator before a new set of pixels
    accumulator.Initialize();

    // Switch between increasing/decreasing scan
    if (m_ScanOrder == INCREASING_ORDER)
      {

      while( !iIt.IsAtEndOfLine() )
        {
	accumulator( iIt.Get() );
	const typename InputImageType::IndexType index = iIt.GetIndex();
	if (outputRegionForThread.IsInside(index))
          {
          outputImage->SetPixel(index, accumulator.GetValue());
          }
	++iIt;
        }
      }
    else
      {
      // DECREASING_ORDER

      iIt.GoToReverseBegin();
      while( !iIt.IsAtReverseEndOfLine() )
        {
        accumulator( iIt.Get() );
	const typename InputImageType::IndexType index = iIt.GetIndex();
	if (outputRegionForThread.IsInside(index))
          {
          outputImage->SetPixel(index, accumulator.GetValue());
          }
	--iIt;
        }
      }

    // Report that the line is finished.
    progress.CompletedPixel();

    // Continue with the next one.
    iIt.NextLine();
  }
}


//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage, class TAccumulator>
TAccumulator
ScanImageFilter<TInputImage,TOutputImage,TAccumulator>
::NewAccumulator(unsigned long size) const
{
  return TAccumulator(size);
}


//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage, class TAccumulator>
void
ScanImageFilter<TInputImage,TOutputImage,TAccumulator>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "ScanDimension: " << m_ScanDimension << std::endl;
  const char *scanOrderString = (m_ScanOrder == INCREASING_ORDER) ?
    "INCREASING_ORDER" : "DECREASING_ORDER";
  os << indent << "ScanOrder: " << scanOrderString << std::endl;
}


} // end namespace itk

#endif // __itkScanImageFilter_cxx
