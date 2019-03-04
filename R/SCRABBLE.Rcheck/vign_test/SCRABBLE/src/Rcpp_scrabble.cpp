// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
using namespace Rcpp;
using namespace Eigen;
using Eigen::SparseMatrix;
using Eigen::MappedSparseMatrix;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::as;
using Eigen::ConjugateGradient;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::Map<Eigen::MatrixXd> MapMatd; // Input: must be double


// [[Rcpp::export]]
SEXP asdgCMatrix_( SEXP XX_ ){

  MapMatd X(as<MapMatd>(XX_));
  SparseMatrix<double> Xsparse = X.sparseView();              // Output: sparse matrix
  S4 Xout(wrap(Xsparse));                                     // Output: as S4 object
  NumericMatrix Xin(XX_);                                     // Copy dimnames
  Xout.slot("Dimnames") = clone(List(Xin.attr("dimnames")));
  return(Xout);
}

// [[Rcpp::export]]
Eigen::VectorXd mSparse(SEXP As, SEXP bs) {
  const MappedSparseMatrix<double> A(as<MappedSparseMatrix<double> >(As));
  const Map<VectorXd> b(as<Map<VectorXd> > (bs));

  VectorXd x=A*b;
  return x;
}

// [[Rcpp::export]]
Eigen::VectorXd mMatrix(SEXP As, SEXP bs) {
  const Map<MatrixXd> A(as<Map<MatrixXd> >(As));
  const Map<VectorXd> b(as<Map<VectorXd> > (bs));

  VectorXd x=A*b;
  return x;
}
// [[Rcpp::export]]
Eigen::VectorXd mSparseT(SEXP As, SEXP bs) {
  const MappedSparseMatrix<double> A(as<MappedSparseMatrix<double> >(As));
  const Map<VectorXd> b(as<Map<VectorXd> > (bs));

  VectorXd x=A.transpose()*b;
  return x;
}
// [[Rcpp::export]]
Eigen::VectorXd mMatrixT(SEXP As, SEXP bs) {
  const Map<MatrixXd> A(as<Map<MatrixXd> >(As));
  const Map<VectorXd> b(as<Map<VectorXd> > (bs));

  VectorXd x=A.transpose()*b;
  return x;
}
// [[Rcpp::export]]
Eigen::MatrixXd getZ(SEXP As, SEXP bs) {
  const Map<MatrixXd> A(as<Map<MatrixXd> >(As));
  double b = Rcpp::as<double>(bs);

  MatrixXd x=b*A.transpose();
  return x;
}

// [[Rcpp::export]]
Eigen::MatrixXd getA(SEXP As, SEXP betas, SEXP gammas, SEXP ns) {
  const Map<MatrixXd> A(as<Map<MatrixXd> >(As));
  double beta = Rcpp::as<double>(betas);
  double gamma = Rcpp::as<double>(gammas);
  double n = Rcpp::as<double>(ns);

  MatrixXd tmp = MatrixXd::Identity(n, n);
  MatrixXd x=beta*(A.transpose()*A)+gamma*tmp;

  return x;
}

// [[Rcpp::export]]
Eigen::MatrixXd getB(SEXP Ds, SEXP Zs, SEXP Ys, SEXP betas) {
  const Map<MatrixXd> D(as<Map<MatrixXd> >(Ds));
  const Map<MatrixXd> Z(as<Map<MatrixXd> >(Zs));
  const MappedSparseMatrix<double> Y(as<MappedSparseMatrix<double> >(Ys));
  double beta = Rcpp::as<double>(betas);

  MatrixXd x=beta*(D.transpose()*Z)+Y;

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd ToDense(SEXP newXs) {
  const MappedSparseMatrix<double> newX(as<MappedSparseMatrix<double> >(newXs));
  return MatrixXd(newX);;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> ToSparse(SEXP newXs) {
  const Map<MatrixXd> newX(as<Map<MatrixXd> >(newXs));
  return newX.sparseView();
}

// [[Rcpp::export]]
Eigen::MatrixXd cDescent(SEXP gamma_Y_B_Lambdas, SEXP As, SEXP zoness, SEXP newXs, SEXP nIters, SEXP error_thresholds){
  double error_threshold = Rcpp::as<double>(error_thresholds);
  int nIter = Rcpp::as<int>(nIters);
  const Map<MatrixXd> gamma_Y_B_Lambda(as<Map<MatrixXd> >(gamma_Y_B_Lambdas));
  const Map<MatrixXd> A(as<Map<MatrixXd> >(As));
  const MappedSparseMatrix<double> zones(as<MappedSparseMatrix<double> >(zoness));
  const Map<MatrixXd> newX(as<Map<MatrixXd> >(newXs));

  int m1 = gamma_Y_B_Lambda.rows();
  int n1 = gamma_Y_B_Lambda.cols();

  int i,j,l;
  double tmp,error;
  /* initalize the outpouts */
  MatrixXd X(m1,n1);
  MatrixXd X1(m1,n1);

  X = newX;
  X1 = newX;
  error = 1;
  l = 1;
  /* Outputs updated */
  while( (l < nIter) && (error > error_threshold)){
    for (i = 0; i < m1; i++){
      for(j = 0; j < n1; j++) {
        tmp = A.row(i)*X.col(j);
        tmp = (gamma_Y_B_Lambda(i,j) - tmp + A(i,i)*X(i,j))/(A(i,i) + zones.coeff(i,j));
        X(i,j) = std::max(tmp,0.);
      }
    }
    error = (X - X1).norm()/(X1.norm());
    X1 = X;
    l = l + 1;
  }
  return X;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> getS(SEXP newXs, SEXP Lambdas, SEXP gammas) {
  const MappedSparseMatrix<double> newX(as<MappedSparseMatrix<double> >(newXs));
  const Map<MatrixXd> Lambda(as<Map<MatrixXd> >(Lambdas));
  double gamma = Rcpp::as<double>(gammas);

  SparseMatrix<double> S=newX+Lambda/gamma;

  return S;
}

// [[Rcpp::export]]
Eigen::MatrixXd getY(SEXP ss, SEXP Us, SEXP Vs) {
  const Map<VectorXd> s(as<Map<VectorXd> >(ss));
  const Map<MatrixXd> U(as<Map<MatrixXd> >(Us));
  const Map<MatrixXd> V(as<Map<MatrixXd> >(Vs));

  MatrixXd S(s.asDiagonal());

  return U*S*(V.transpose());
}

// [[Rcpp::export]]
double calculateError(SEXP Xs, SEXP newXs, SEXP m1s, SEXP n1s) {
  const MappedSparseMatrix<double> newX(as<MappedSparseMatrix<double> >(newXs));
  const MappedSparseMatrix<double> X(as<MappedSparseMatrix<double> >(Xs));
  int m1 = Rcpp::as<int>(m1s);
  int n1 = Rcpp::as<int>(n1s);
  const MatrixXd C = MatrixXd::Ones(m1,n1);

  MatrixXd tmp1 = MatrixXd(newX) + C;
  MatrixXd tmp2 = MatrixXd(X)+C;

  double error = ((tmp1.log() - tmp2.log())/log(10.0)).norm()/(m1*n1);

  return error;
}

// [[Rcpp::export]]
Eigen::MatrixXd updateLambda(SEXP Lambdas, SEXP newXs, SEXP newYs, SEXP gammas) {
  const MappedSparseMatrix<double> newX(as<MappedSparseMatrix<double> >(newXs));
  const Map<MatrixXd> Lambda(as<Map<MatrixXd> >(Lambdas));
  const Map<MatrixXd> newY(as<Map<MatrixXd> >(newYs));
  double gamma = Rcpp::as<double>(gammas);

  MatrixXd X=Lambda + gamma*(newX - newY);

  return X;
}

// [[Rcpp::export]]
Eigen::MatrixXd recoverData(SEXP newXs) {
  const MappedSparseMatrix<double> newX(as<MappedSparseMatrix<double> >(newXs));

  MatrixXd X=(MatrixXd(newX)).transpose();

  return X;
}
