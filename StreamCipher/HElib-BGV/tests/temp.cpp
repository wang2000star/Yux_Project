uint64_t nonce = 123456789;
  size_t size = ciphertexts.size();
  std::cout << "size: " << size << std::endl;
  std::cout << "cipher_size: " << params.cipher_size << std::endl;
  size_t num_block = ceil((double)size / params.cipher_size);

  Pasta2 pasta2(instance);

  std::vector<std::vector<uint64_t>> fixed_mat1 = instance.matrix;
  std::vector<std::vector<uint64_t>> fixed_mat2 = instance.matrix;

  const helib::PAlgebra& zMStar = context->getZMStar();
  long root;
  long Para_m = 65536;
  if (Para_m == 65536) {
    root = 3;
  }
  if (Para_m == 32768) {
    root = 9;
  }
  if (Para_m == 16384) {
    root = 81;
  }
  long Para_p = 65537;
  long slots_len = 32768;
  long plain_size = 128;
  long key_size = 256;
  long plain_size_square = plain_size * plain_size;
  // 初始化Cmodulus对象
  helib::Cmodulus cmodulus(zMStar, Para_p, root);
  zz_pX encodedtemp;
  // 加密对称密钥
  Ctxt tmpctxt(*he_pk);
  std::vector<Ctxt> encrypedSymKey(key_size, tmpctxt);
  std::vector<Ctxt> encrypedKeyStream(key_size, tmpctxt);
  vec_long slotsData;
  slotsData.SetLength(slots_len);
  ZZX encodeddata;
  for (int i = 0; i < key_size; i++) {
    for (int j = 0; j < slots_len; j++) {
      slotsData[j] = secret_key[i];
    }
    cmodulus.iFFT(encodedtemp, slotsData);
    convert(encodeddata, encodedtemp);
    he_pk->Encrypt(encrypedKeyStream[i], encodeddata);
  }
  std::cout << "KeyStream 加密完成！" << std::endl;
  print_noise(encrypedKeyStream);
  std::vector<Ctxt> temp = encrypedKeyStream;

  auto start_time = std::chrono::high_resolution_clock::now();
  auto start_server_off = std::chrono::high_resolution_clock::now();
  // num_block=slots_len
  //  XOF
  std::vector<vec_long> xof_mat1(plain_size_square);
  for (int i = 0; i < plain_size_square; i++) {
    xof_mat1[i].SetLength(slots_len);
  }
  ///(128*128, std::vector<long>(slots_len, 0));
  std::vector<vec_long> xof_mat2 = xof_mat1;
  std::vector<vec_long> xof_rc(key_size);
  for (int i = 0; i < key_size; i++) {
    xof_rc[i].SetLength(slots_len);
  }
  auto start_xof = std::chrono::high_resolution_clock::now();
  for (uint64_t b = 0; b < num_block; b++) {
    pasta2.init_shake(nonce, b);
    auto mat1 = pasta2.get_random_matrixL();
    auto mat2 = pasta2.get_random_matrixR();
    auto rc = pasta2.get_rc_vec(2 * PASTA2_T);
    for (int i = 0; i < plain_size; i++) {
      for (int j = 0; j < plain_size; j++) {
        xof_mat1[i * plain_size + j][b] = mat1[i][j];
        xof_mat2[i * plain_size + j][b] = mat2[i][j];
      }
      xof_rc[i][b] = rc[i];
      xof_rc[i + plain_size][b] = rc[i + plain_size];
    }
  }
  auto end_xof = std::chrono::high_resolution_clock::now();
  double xof_time = std::chrono::duration_cast<std::chrono::duration<double>>(
                        end_xof - start_xof)
                        .count();
  std::cout << "XOF生成时间：" << xof_time << std::endl;
  // 对
  std::vector<zzX> encodedxof_mat1(plain_size_square);
  std::vector<zzX> encodedxof_mat2(plain_size_square);
  std::vector<ZZX> encodedxof_rc(key_size);
  auto start_encode = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < plain_size_square; i++) {
    cmodulus.iFFT(encodedtemp, xof_mat1[i]);
    convert(encodedxof_mat1[i], encodedtemp);
    cmodulus.iFFT(encodedtemp, xof_mat2[i]);
    convert(encodedxof_mat1[i], encodedtemp);
  }
  for (int i = 0; i < key_size; i++) {
    cmodulus.iFFT(encodedtemp, xof_rc[i]);
    convert(encodedxof_rc[i], encodedtemp);
  }
  auto end_encode = std::chrono::high_resolution_clock::now();
  double encode_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(end_encode -
                                                                start_encode)
          .count();
  std::cout << "编码时间：" << encode_time << std::endl;

  // 下面执行随机线性层
  Ctxt temp1 = tmpctxt;
  auto start_random = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < plain_size; i++) {
    encrypedKeyStream[i] = temp[0];
    encrypedKeyStream[i].multByConstant(encodedxof_mat1[i * plain_size + 0]);
    for (int j = 1; j < plain_size; j++) {
      temp1 = temp[j];
      temp1.multByConstant(encodedxof_mat1[i * plain_size + j]);
      encrypedKeyStream[i] += temp1;
    }
    encrypedKeyStream[i].addConstant(encodedxof_rc[i]);
  }
  for (int i = 0; i < plain_size; i++) {
    encrypedKeyStream[i + plain_size] = temp[plain_size];
    encrypedKeyStream[i + plain_size].multByConstant(
        encodedxof_mat1[i * plain_size + 0]);
    for (int j = 1; j < plain_size; j++) {
      temp1 = temp[plain_size + j];
      temp1.multByConstant(encodedxof_mat1[i * plain_size + j]);
      encrypedKeyStream[i] += temp1;
    }
    encrypedKeyStream[i + plain_size].addConstant(
        encodedxof_rc[i + plain_size]);
  }
  temp = encrypedKeyStream;
  for (int i = 0; i < key_size; i++) {
    encrypedKeyStream[i] += temp[i];
    encrypedKeyStream[i] += temp[(i + plain_size) % key_size];
  }
  auto end_random = std::chrono::high_resolution_clock::now();
  double random_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(end_random -
                                                                start_random)
          .count();
  // 随机层结束
  print_noise(encrypedKeyStream);
  auto start_round = std::chrono::high_resolution_clock::now();
  auto end_round = std::chrono::high_resolution_clock::now();
  double round_time = 0;
  for (int r = 0; r < 3; r++) {
    std::cout << "round: " << r << std::endl;
    start_round = std::chrono::high_resolution_clock::now();
    if (r == 2)
      for (int i = 0; i < key_size; i++) {
        encrypedKeyStream[i].cube();
      }
    else {
      temp = encrypedKeyStream;
      for (int i = 1; i < plain_size; i++) {
        temp[i - 1].square();
        encrypedKeyStream[i] += temp[i - 1];
        temp[i + plain_size - 1].square();
        encrypedKeyStream[i + plain_size] += temp[i + plain_size - 1];
      }
    }
    temp = encrypedKeyStream;
    for (int i = 0; i < plain_size; i++) {
      encrypedKeyStream[i] = temp[0];
      encrypedKeyStream[i].multByConstant(static_cast<ZZX>(fixed_mat1[i][0]));
      for (int j = 1; j < plain_size; j++) {
        temp1 = temp[j];
        temp1.multByConstant(static_cast<ZZX>(fixed_mat1[i][j]));
        encrypedKeyStream[i] += temp1;
      }
    }
    for (int i = 0; i < plain_size; i++) {
      encrypedKeyStream[i + plain_size] = temp[plain_size];
      encrypedKeyStream[i + plain_size].multByConstant(
          static_cast<ZZX>(fixed_mat2[i][0]));
      for (int j = 1; j < plain_size; j++) {
        temp1 = temp[plain_size + j];
        temp1.multByConstant(static_cast<ZZX>(fixed_mat2[i][j]));
        encrypedKeyStream[i + plain_size] += temp1;
      }
    }
    for (int i = 0; i < key_size; i++) {
      encrypedKeyStream[i].addConstant(to_ZZX(r));
    }
    temp = encrypedKeyStream;
    for (int i = 0; i < key_size; i++) {
      encrypedKeyStream[i] += temp[i];
      encrypedKeyStream[i] += temp[(i + plain_size) % key_size];
    }
    end_round = std::chrono::high_resolution_clock::now();
    round_time += std::chrono::duration_cast<std::chrono::duration<double>>(
                      end_round - start_round)
                      .count();
    print_noise(encrypedKeyStream);
  }
  auto end_server_off = std::chrono::high_resolution_clock::now();
  double server_off_time = xof_time + encode_time + random_time + round_time;
  std::vector<Ctxt> TruncencryptedKeyStream(
      encrypedKeyStream.begin(), encrypedKeyStream.begin() + plain_size);
  std::vector<Ctxt> encrypedPlainStream = TruncencryptedKeyStream;
  print_noise(encrypedPlainStream);
  std::vector<vec_long> CipherStream(plain_size);
  for (int i = 0; i < plain_size; i++) {
    CipherStream[i].SetLength(slots_len);
  }
  random_device rd;
  for (int i = 0; i < plain_size; i++) {
    for (int j = 0; j < slots_len; j++) {
      CipherStream[i][j] = rd() % 65537;
    }
  }
  // 对CipherStream进行编码
  std::vector<ZZX> encodedCipherStream(plain_size);
  auto start_ServerOn = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < plain_size; i++) {
    cmodulus.iFFT(encodedtemp, CipherStream[i]);
    convert(encodedCipherStream[i], encodedtemp);
  }
  for (int i = 0; i < plain_size; i++) {
    encrypedPlainStream[i].negate();
    encrypedPlainStream[i].addConstant(encodedCipherStream[i]);
  }
  auto end_ServerOn = std::chrono::high_resolution_clock::now();
  double ServerOn_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(end_ServerOn -
                                                                start_ServerOn)
          .count();
  print_noise(encrypedPlainStream);
  auto end_time = std::chrono::high_resolution_clock::now();
  std::cout << "ServerOff time: " << server_off_time << std::endl;
  std::cout << "ServerOn time: " << ServerOn_time << std::endl;
  std::cout << "ServerTotal time: "
            << std::chrono::duration_cast<std::chrono::duration<double>>(
                   end_time - start_time)
                   .count()
            << std::endl