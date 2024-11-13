#include <helib/helib.h>
#include <vector>
#include <iostream>

using namespace std;
using namespace helib;

// Function to create a MatMul1D object from a given matrix
class MyMatMul1D : public MatMul1D_derived<PA_GF2>
{
private:
  const EncryptedArray& ea;
  vector<vector<long>> matrix;

public:
  MyMatMul1D(const EncryptedArray& _ea, const vector<vector<long>>& _matrix) : ea(_ea), matrix(_matrix) {}

  virtual const EncryptedArray& getEA() const override { return ea; }

  virtual long getDim() const override { return 0; }

  virtual bool multipleTransforms() const override { return false; }

  virtual bool get(RX& out, long i, long j, long k) const override
  {
    if (i < matrix.size() && j < matrix[i].size()) {
      out = matrix[i][j];
      return false;
    }
    return true;
  }
};

int main()
{
  // Initialize context
  long p = 2;    // Plaintext prime modulus
  long r = 1;    // Lifting
  long L = 16;   // Number of bits of the modulus chain
  long c = 2;    // Columns in the key-switching matrix
  long k = 80;   // Security parameter
  long s = 0;    // Minimum number of slots
  long d = 1;    // Degree of the field extension

  Context context = ContextBuilder<BGV>()
                        .m(8192)
                        .p(p)
                        .r(r)
                        .bits(L)
                        .c(c)
                        .build();

  // Generate secret key
  SecKey secretKey(context);
  secretKey.GenSecKey();

  // Generate public key
  const PubKey& publicKey = secretKey;

  // Get the EncryptedArray of the context
  const EncryptedArray& ea = context.getEA();

  // Define matrix M and vector v
  vector<vector<long>> M = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9}};
  vector<long> v = {1, 2, 3};

  // Encrypt the vector
  PtxtArray ptxt_v(ea);
  ptxt_v.encode(v);
  Ctxt ctxt_v(publicKey);
  ptxt_v.encrypt(ctxt_v);

  // Create MatMul1D object
  MyMatMul1D matMul(ea, M);

  // Perform matrix-vector multiplication
  ctxt_v *= matMul;

  // Decrypt the result
  PtxtArray ptxt_result(ea);
  ptxt_result.decrypt(ctxt_v, secretKey);

  // Print the result
  vector<long> result;
  ptxt_result.decode(result);
  cout << "Result: ";
  for (const auto& val : result) {
    cout << val << " ";
  }
  cout << endl;

  return 0;
}
