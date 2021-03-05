import 'dart:async';
import 'dart:io';

import '../lib/lzma.dart';

Future main(List<String> args) async {
  if (args.length != 3) {
    print('Usage: compress input output');
    return;
  }

  final inFile = File(args[0]);
  final outFile = File(args[1]);
  final decodeOutFile = File(args[2]);

  final List<int?> encoded = lzma.encode(await inFile.readAsBytes());
  await outFile.writeAsBytes(encoded.map((e) => e!).toList());

  final List<int?> decode = lzma.decode(await outFile.readAsBytes());
  await decodeOutFile.writeAsBytes(decode.map((e) => e!).toList());
}
