/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 2.0.4
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.sbml.libsbml;

public class SBMLNamespacesList {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  public SBMLNamespacesList(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  public static long getCPtr(SBMLNamespacesList obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        libsbmlJNI.delete_SBMLNamespacesList(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public SBMLNamespacesList() {
    this(libsbmlJNI.new_SBMLNamespacesList(), true);
  }

  public void add(SBMLNamespaces item) {
    libsbmlJNI.SBMLNamespacesList_add(swigCPtr, this, SBMLNamespaces.getCPtr(item), item);
  }

  public SBMLNamespaces get(long n) {
    long cPtr = libsbmlJNI.SBMLNamespacesList_get(swigCPtr, this, n);
    return (cPtr == 0) ? null : new SBMLNamespaces(cPtr, false);
  }

  public void prepend(SBMLNamespaces item) {
    libsbmlJNI.SBMLNamespacesList_prepend(swigCPtr, this, SBMLNamespaces.getCPtr(item), item);
  }

  public SBMLNamespaces remove(long n) {
    long cPtr = libsbmlJNI.SBMLNamespacesList_remove(swigCPtr, this, n);
    return (cPtr == 0) ? null : new SBMLNamespaces(cPtr, false);
  }

  public long getSize() {
    return libsbmlJNI.SBMLNamespacesList_getSize(swigCPtr, this);
  }

}
