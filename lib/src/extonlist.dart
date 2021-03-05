extension ListputIf<E> on List<E> {
  void replaceOrAdd(E val, int index) {
    if (this.length > index) {
      print("this.length > index $index");
      this[index] = val;
    } else if (length == index) {
      print("length==index $index");

      this.add(val);
    } else {
      print("else ");

      // non supported right now
    }
  }

  E? elementAtOrNull(int index){
      if(index>=this.length) return null;
      return this[index];
  }
}
