These are some sample scripts for doing volume rendering with yt.
These were originally written with the WD merger problem in mind.

  -- vol-wd.py :

     this script does basic volume rendering using the perspective
     lens.

  -- vol-wd-spherical.py :

     this script does a spherical projection (all 4 pi steradians)
     with the camera set at the origin.  This can be used to create
     360 degree videos for uploading to youtube.

     Note to prepare the video for youtube, follow the instructions
     here:

     https://support.google.com/youtube/answer/6178631?hl=en

     In particular, you need to add metadata to the video, which can be
     done with the google spatial media metadata injector:

     https://github.com/google/spatial-media/blob/master/spatialmedia/README.md


