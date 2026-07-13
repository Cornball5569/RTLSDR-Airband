                    // If squelch is open / opening and using I/Q, then cleanup the signal and possibly update squelch.
                     if (fparms->squelch.should_filter_sample() && channel->needs_raw_iq) {
                         // remove phase rotation introduced by FFT sliding window
                         float swf, cwf, re_tmp, im_tmp;
                         sincosf_lut(channel->dm_phi, &swf, &cwf);
                         multiply(real, imag, cwf, -swf, &re_tmp, &im_tmp);
                         channel->dm_phi += channel->dm_dphi;
                         channel->dm_phi &= 0xffffff;

                         // apply lowpass filter, will be a no-op if not configured
                         fparms->lowpass_filter.apply(re_tmp, im_tmp);

                         // apply bandpass filter for audio, will be a no-op if not configured
                         fparms->bandpass_filter.apply(re_tmp, im_tmp);

                         // update I/Q and wave
                         real = re_tmp;
                         imag = im_tmp;
                         channel->wavein[j] = sqrt(real * real + imag * imag);

                         // update squelch post-cleanup
                         if (fparms->lowpass_filter.enabled() || fparms->bandpass_filter.enabled()) {
                             fparms->squelch.process_filtered_sample(channel->wavein[j]);
                         }
                     }
