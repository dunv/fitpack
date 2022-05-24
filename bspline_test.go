package fitpack

import (
	"testing"

	"github.com/stretchr/testify/require"
)

var points = [][]float64{
	{1.889, 7.100},
	{4.164, 7.100},
	{5.743, 7.100},
	{8.439, 7.100},
	{9.876, 7.100},
	{11.000, 7.000},
	{12.200, 6.313},
	{12.750, 4.444},
	{12.500, 3.000},
	{11.000, 1.900},
	{9.793, 1.800},
	{6.906, 1.800},
	{4.836, 1.800},
	{3.530, 1.800},
	{1.553, 1.800},
	{0.100, 2.094},
	{-0.900, 3.446},
	{-1.100, 4.115},
	{-1.000, 5.221},
	{-0.250, 6.586},
	{0.800, 7.000},
	{1.889, 7.100},
}

func TestSpline(t *testing.T) {
	spline, err := NewBSpline(points, 0.2)
	require.NoError(t, err)

	pt, err := spline.Splev([]float64{.3})
	require.NoError(t, err)

	require.Len(t, pt, 1)
	require.Len(t, pt[0], 2)
	require.Equal(t, 11.946815183056724, pt[0][0])
	require.Equal(t, 6.417278208932581, pt[0][1])

}
