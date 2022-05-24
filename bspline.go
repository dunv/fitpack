package fitpack

// #cgo LDFLAGS: -L/usr/bin/gcc -lgfortran -Wl,--allow-multiple-definition -lm
// #cgo CFLAGS: -w
// void parcur_( int*, int*, int*, int*, double*, int*, double*, double*, double*, double*, int*, double*, int*, int*, double*, int*, double*, double*, double*, int*, int*, int* );
// void splev_( double*, int*, double*, int*, double*, double*, int*, int*, int* );
// void splder_( double*, int*, double*, int*, int*, double*, double*, int*, int*, double*, int* );
import "C"

import (
	"errors"
	"fmt"
)

type BSpline struct {
	t []float64
	c [][]float64
	k int32
	u []float64
}

func NewBSpline(points [][]float64, s float64) (*BSpline, error) {
	iopt := int32(0)
	ipar := int32(0)
	idim := int32(len(points[0]))
	m := int32(len(points))
	u := make([]float64, m)

	x := []float64{}
	for _, coord := range points {
		x = append(x, coord[0], coord[1])
	}

	mx := int32(len(x))
	w := make([]float64, 0)
	for i := int32(0); i < m; i++ {
		w = append(w, 1)
	}

	ub := 0.0
	ue := 1.0
	k := int32(3)
	nest := m + 2*k // always enough according to docu m+k+1
	var n int32
	t := make([]float64, nest)
	nc := nest * idim
	c := make([]float64, nc)
	var fp float64
	lwrk := m*(k+1) + nest*(6+idim+3*k)
	wrk := make([]float64, lwrk)

	iwrk := make([]int32, nest)
	var ier int32

	if !(-1 <= iopt) || !(iopt <= 1) {
		return nil, errors.New("err -1<=iopt<=1 does not hold")
	}
	if !(0 <= k) || !(k <= 5) {
		return nil, errors.New("err 0<=k<=5 does not hold")
	}
	if !(m > k) {
		return nil, errors.New("err m>k does not hold")
	}
	if !(nest > 2*k+2) {
		return nil, errors.New("err nest > 2*k+2 does not hold")
	}
	if !(int32(len(w)) == m) {
		return nil, errors.New("err len(w) == m does not hold")
	}
	for _, val := range w {
		if !(val > 0) {
			return nil, errors.New("err w(i)>0,i=1,2,...,m does not hold")
		}
	}
	if !(0 <= ipar) || !(ipar <= 1) {
		return nil, errors.New("err 0<=ipar<=1 does not hold")
	}
	if !(0 <= idim) || !(idim <= 10) {
		return nil, errors.New("err 0<=ipar<=10 does not hold")
	}
	if !(lwrk >= (k+1)*m+nest*(6+idim+3*k)) {
		return nil, errors.New("err lwrk >= (k+1)*m+nest*(6+idim+3*k) does not hold")
	}
	if !(nc >= nest*idim) {
		return nil, errors.New("err nc >= nest*idim does not hold")
	}

	C.parcur_(
		(*C.int)(&iopt),
		(*C.int)(&ipar),
		(*C.int)(&idim),
		(*C.int)(&m),
		(*C.double)(&u[0]),
		(*C.int)(&mx),
		(*C.double)(&x[0]),
		(*C.double)(&w[0]),
		(*C.double)(&ub),
		(*C.double)(&ue),
		(*C.int)(&k),
		(*C.double)(&s),
		(*C.int)(&nest),
		(*C.int)(&n),
		(*C.double)(&t[0]),
		(*C.int)(&nc),
		(*C.double)(&c[0]),
		(*C.double)(&fp),
		(*C.double)(&wrk[0]),
		(*C.int)(&lwrk),
		(*C.int)(&iwrk[0]),
		(*C.int)(&ier),
	)
	if ier != 0 {
		return nil, fmt.Errorf("err fitpack/parcur returns non-zero exitCode (%d)", ier)
	}

	t = t[:n]
	cShaped := [][]float64{}
	for j := int32(0); j < idim; j++ {
		cShaped = append(cShaped, make([]float64, n-k-1))
		for i := int32(0); i < n-k-1; i++ {
			cShaped[j][i] = c[n*(j)+i]
		}
	}

	spline := &BSpline{
		t: t,
		c: cShaped,
		k: k,
		u: u,
	}

	return spline, nil
}

// void splev_( double*, int*,  double*, int*, double*, double*, int*,  int*, int* );
// 				t		 len(t) c		 k	   x		y		 len(x) e	  ier
func (s *BSpline) Splev(pos []float64) ([][]float64, error) {
	n := int32(len(s.t))
	if len(pos) < 1 {
		return nil, errors.New("at least one point is required for calculation")
	}
	m := int32(len(pos))
	e := int32(0)

	output := []float64{}
	for index := range s.c {
		y := make([]float64, m)
		ier := int32(0)
		C.splev_(
			(*C.double)(&s.t[0]),
			(*C.int)(&n),
			(*C.double)(&s.c[index][0]),
			(*C.int)(&s.k),
			(*C.double)(&pos[0]),
			(*C.double)(&y[0]),
			(*C.int)(&m),
			(*C.int)(&e),
			(*C.int)(&ier),
		)
		if ier != 0 {
			return nil, errors.New("err fitpack/splev returns non-zero exitCode")
		}
		output = append(output, y...)
	}

	shapedOutput := [][]float64{}
	for i := 0; i < len(output)/2; i++ {
		shapedOutput = append(shapedOutput, []float64{output[i], output[i+len(output)/2]})
	}

	return shapedOutput, nil
}

// void splder_( double*, int*, double*, int*, int*, double*, double*, int*, int*, double*, int* );
// 				t		n		c		k		nu	 x			y	  len(x) e     wrk      ier
func (s *BSpline) Splder(pos []float64, der int32) ([][]float64, error) {
	n := int32(len(s.t))
	m := int32(len(pos))
	e := int32(0)
	nu := int32(1)

	output := []float64{}
	for index := range s.c {
		y := make([]float64, m)
		wrk := make([]float64, n)
		ier := int32(0)
		C.splder_(
			(*C.double)(&s.t[0]),
			(*C.int)(&n),
			(*C.double)(&s.c[index][0]),
			(*C.int)(&s.k),
			(*C.int)(&nu),
			(*C.double)(&pos[0]),
			(*C.double)(&y[0]),
			(*C.int)(&m),
			(*C.int)(&e),
			(*C.double)(&wrk[0]),
			(*C.int)(&ier),
		)
		if ier != 0 {
			return nil, errors.New("err fitpack/spder returns non-zero exitCode")
		}
		output = append(output, y...)
	}
	shapedOutput := [][]float64{}
	for i := 0; i < len(output)/2; i++ {
		shapedOutput = append(shapedOutput, []float64{output[i], output[i+len(output)/2]})
	}
	return shapedOutput, nil
}
