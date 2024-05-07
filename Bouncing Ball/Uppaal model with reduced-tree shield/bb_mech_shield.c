int shield(double E, double v, double p_thresh) {
   if (E <= 48.0) {
    if (E <= 16.0) {
      if (E <= 8.0) {
        if (E <= 4.0) {
          return 0;
        } else {
          if (v <= -4.0) {
            return 3;
          } else {
            if (v <= 4.0) {
              if (p_thresh <= 0.5) {
                return 0;
              } else {
                return 3;
              }
            } else {
              return 3;
            }
          }
        }
      } else {
        if (v <= -6.0) {
          return 3;
        } else {
          if (v <= 6.0) {
            if (p_thresh <= 0.5) {
              return 0;
            } else {
              return 3;
            }
          } else {
            return 3;
          }
        }
      }
    } else {
      if (v <= 2.0) {
        if (v <= -8.0) {
          if (v <= -10.0) {
            return 3;
          } else {
            if (E <= 32.0) {
              return 3;
            } else {
              if (p_thresh <= 0.5) {
                return 0;
              } else {
                return 3;
              }
            }
          }
        } else {
          if (p_thresh <= 0.5) {
            if (E <= 40.0) {
              return 0;
            } else {
              if (v <= -2.0) {
                return 0;
              } else {
                return 3;
              }
            }
          } else {
            if (E <= 36.0) {
              return 3;
            } else {
              if (v <= -2.0) {
                if (E <= 40.0) {
                  return 3;
                } else {
                  if (v <= -4.0) {
                    return 3;
                  } else {
                    return 0;
                  }
                }
              } else {
                return 0;
              }
            }
          }
        }
      } else {
        if (p_thresh <= 0.5) {
          if (E <= 40.0) {
            if (v <= 8.0) {
              return 0;
            } else {
              if (E <= 32.0) {
                return 3;
              } else {
                if (v <= 10.0) {
                  return 0;
                } else {
                  return 3;
                }
              }
            }
          } else {
            if (v <= 10.0) {
              return 0;
            } else {
              return 3;
            }
          }
        } else {
          if (E <= 40.0) {
            return 3;
          } else {
            if (v <= 4.0) {
              return 0;
            } else {
              return 3;
            }
          }
        }
      }
    }
  } else {
    if (v <= 2.0) {
      if (p_thresh <= 0.5) {
        if (E <= 56.0) {
          if (v <= -12.0) {
            return 3;
          } else {
            if (v <= -4.0) {
              return 0;
            } else {
              return 3;
            }
          }
        } else {
          if (E <= 72.0) {
            if (v <= -12.0) {
              return 3;
            } else {
              if (v <= -6.0) {
                return 0;
              } else {
                return 3;
              }
            }
          } else {
            return 3;
          }
        }
      } else {
        if (v <= -2.0) {
          if (E <= 56.0) {
            if (v <= -6.0) {
              return 3;
            } else {
              return 0;
            }
          } else {
            if (v <= -8.0) {
              return 3;
            } else {
              if (E <= 72.0) {
                return 0;
              } else {
                return 3;
              }
            }
          }
        } else {
          if (E <= 68.0) {
            if (E <= 64.0) {
              return 0;
            } else {
              if (v <= 0.0) {
                return 0;
              } else {
                return 1;
              }
            }
          } else {
            if (E <= 72.0) {
              if (v <= 0.0) {
                return 1;
              } else {
                return 3;
              }
            } else {
              return 3;
            }
          }
        }
      }
    } else {
      if (v <= 4.0) {
        if (p_thresh <= 0.5) {
          return 3;
        } else {
          if (E <= 64.0) {
            if (E <= 60.0) {
              return 0;
            } else {
              return 1;
            }
          } else {
            return 3;
          }
        }
      } else {
        if (E <= 52.0) {
          if (v <= 6.0) {
            return 0;
          } else {
            if (v <= 12.0) {
              if (p_thresh <= 0.5) {
                return 0;
              } else {
                return 3;
              }
            } else {
              return 3;
            }
          }
        } else {
          if (E <= 64.0) {
            if (v <= 6.0) {
              if (p_thresh <= 0.5) {
                return 3;
              } else {
                if (E <= 60.0) {
                  return 1;
                } else {
                  return 3;
                }
              }
            } else {
              return 3;
            }
          } else {
            return 3;
          }
        }
      }
    }
  }

}